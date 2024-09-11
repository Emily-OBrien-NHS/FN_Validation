from data_cleaning_and_transformation import sort_events
from pathlib import Path
from pm4py.objects.conversion.log import converter as log_converter
from pm4py.algo.discovery.dfg import algorithm as dfg_discovery
from pm4py.visualization.dfg import visualizer as dfg_visualizer
from graphviz import Digraph
import pandas as pd
import numpy as np

################################################################################
#------------------pathway definitions supporting functions--------------------#
################################################################################

def add_reset_transitions(events_data):
    #Function to add the count, total and percentage columns to the events
    #data frame to work out pathways.  If already done, removes the original
    #columns
    transitions = events_data.copy()
    transitions = transitions.drop(["Count", "Total", "Percentage"], axis=1,
                                   errors="ignore", inplace=False)
    # Calculate Next Event for each Patient
    transitions["Next Event (Pathway)"] = (transitions.groupby("VisitId")
                                           ["Event (Pathway)"].shift(-1))
    # Calculate the count of patients that make each transition
    count_event_pairs = (transitions
                         .groupby(["Event (Pathway)", "Next Event (Pathway)"],
                                  as_index=False).agg({"VisitId": "count"})
                                  .rename(columns = {"VisitId": "Count"})
                                  .drop_duplicates())
    #calculate the total number of patients that have each event
    count_events = (transitions
                    .loc[~transitions["Next Event (Pathway)"].isnull()]
                    .groupby(["Event (Pathway)"], as_index=False)
                    .agg({"VisitId": "count"})
                    .rename(columns={"VisitId": "Total"}).drop_duplicates())
    #merge counts onto events data
    transitions = (transitions
                   .merge(count_event_pairs,
                          on=["Event (Pathway)", "Next Event (Pathway)"],
                          how="left")
                   .merge(count_events, on=["Event (Pathway)"], how="left"))
    #calculate the percentage of patients that make each transition from each
    #event
    transitions["Percentage"] = (100 * transitions["Count"]
                                 / transitions["Total"])
    return transitions

def generate_pathway_definitions(events_data, include_spawn_end_events):
    #Function to create the pathway definition
    #----
    #If spawn and end events were included, remove them here for next steps.
    if include_spawn_end_events:
        events_data = events_data.loc[~events_data["EventName"]
                                      .isin(["Spawn", "Removed"])].copy()
    #Remove duplicates from to events, remove unecessary columns and rename.
    pathway_definition = (events_data
                         .drop_duplicates(["Event (Pathway)",
                                            "Next Event (Pathway)"])
                         .drop(columns=["Count", "Total"])
                         .rename(columns={"Event (Pathway)": "From Process",
                                         "Next Event (Pathway)": "To Process"}))
    #Add empty priority and notes columns
    pathway_definition["(Consequent Priority)"] = None
    pathway_definition["Notes"] = None
    #filter to transitions that have a greater than 0 percentage
    pathway_definition = (pathway_definition
                          .loc[pathway_definition["Percentage"] >= 0,
                               ["From Process", "To Process",
                                "(Consequent Priority)", "Percentage",
                                "Notes"]].copy())
    return pathway_definition

def add_obs_repeat_splits(pathway_definition, obs_splits):
    #Function to add in the required 'kick off' events to the events dataframe
    #which are needed to start repeated obs events
    #----
    #Filter the pathway to the triage events that need the extra event to
    #trigger different obs repetitions. i[0]==i[0] to remove nan in From event
    #where kick off events are not required. Create df of these and remove them
    #from the pathway definitions.
    triage_events_for_obs = set([i[0] for i in obs_splits if i[0]==i[0]])
    triage_events_mask = (pathway_definition['From Process']
                          .isin(triage_events_for_obs))
    triage_events = pathway_definition.loc[triage_events_mask].copy()
    pathway_definition = pathway_definition.loc[~triage_events_mask].copy()
    #Loop through the list of splits, and create a new dataframe with triaged
    #events going to the triaged obs events with their probabilities.
    triage_obs = []
    for triage_values in obs_splits:
        #Only add events that require a kick off event (anywhere that
        #triage_values[0] (From Process) is not nan)
        if triage_values[0]==triage_values[0]:
            triage_obs.append(
                [triage_values[0], triage_values[1], None,
                 float(triage_values[2]),
                'Process to kick off repeated obs for different timings'])
    triage_obs_events = pd.DataFrame(triage_obs, columns=triage_events.columns)
    #Create a dataframe of the inserted triage obs events to the original triage
    #'To Process'.
    new_links = triage_obs_events[['From Process', 'To Process']]
    new_links.columns = ['Original Process', 'From Process']
    triage_events = triage_events.rename(columns=
                                         {'From Process':'Original Process'})
    new_links = (new_links.merge(triage_events, on='Original Process')
                 .drop('Original Process', axis=1))
    #concat the triage to new triage obs processes dataframe and the new triage
    #obs processes to the original 'To Process' onto the pathway definition.
    pathway_definition = pd.concat([triage_obs_events, new_links,
                                    pathway_definition])
    return pathway_definition

def create_process_recurrence(obs_splits):
    #Function to create the process recurrence outputs.
    #----
    #Create the process recurrence triggers output.
    process_recur_triggers = pd.DataFrame([lst[1:] for lst in obs_splits],
                                  columns=['Trigger Process (In Pathway)',
                                          'Probability (%)',
                                          'Recurrent Process (Not In Pathway)'])
    #As we've added in kick-off events in the pathway to do probability splits,
    #all probabilites here are 100%
    process_recur_triggers['Probability (%)'] = 100
    process_recur_triggers['Notes'] = np.nan
    #Rearrange columns
    process_recur_triggers = process_recur_triggers[[
                             'Trigger Process (In Pathway)',
                             'Recurrent Process (Not In Pathway)',
                             'Probability (%)', 'Notes']].copy()
    #Create the process recurrence output.
    process_recur = pd.DataFrame(process_recur_triggers
                                 ['Recurrent Process (Not In Pathway)']
                                 .rename('Recurrent Process'))
    process_recur['Recurrence Mean'] = (process_recur['Recurrent Process'].str
                                        .extract(r'(\d+)').astype(int))
    process_recur['StdDev'] = (round(0.1 * process_recur['Recurrence Mean'])
                               .astype(int))
    process_recur['Min'] = 5
    process_recur['Max'] = np.nan
    process_recur['Notes'] = np.nan
    return process_recur_triggers, process_recur

def output_transition_viz(pathway_definition, filepath, name):
    #Function to create the pathway diagrams as webpage for all pathways
    def get_colour(node):
        #Function for the colour of each node
        colour = ('lightcoral' if ('Arrival' in node) or ('Walk-In' in node)
                  else 'lightgreen'
                  if ('Admitted' in node and not 'Wait' in node)
                  or ('Discharged' in node and not 'Wait' in node)
                  else 'aliceblue')
        return colour
    def get_shape(node):
        #Function for the shape of each node
        shape = ('invhouse' if ('Arrival' in node) or ('Walk-In' in node)
                  else 'house'
                  if ('Admitted' in node and not 'Wait' in node)
                  or ('Discharged' in node and not 'Wait' in node)
                  else 'oval')
        return shape
    #Sort the pathway definition
    pathway_definition = (pathway_definition
                          .sort_values(["To Process", "Percentage"]))
    #Start setting up the plot.
    graph = Digraph(name, engine='dot')
    graph.attr(fontsize='16', fontname = "Segoe UI", rankdir='LR')
    #get list of all the nodes.
    nodes = set(pd.concat([pathway_definition["To Process"],
                           pathway_definition["From Process"]]).to_list())
    #for each node, get the colour and shape
    for node in nodes:
        graph.node(node, style='filled', group='none',
                   fillcolor=get_colour(node), shape=get_shape(node),
                   fontsize='22', penwidth='2')
    #plot each arrow with the percentage as text, and the width based on the
    #proportion of patients going down that path.
    for _, row in pathway_definition.iterrows():
        graph.edge(row["From Process"], row["To Process"],
                   "{:.1f}%".format(row["Percentage"]),
                   penwidth = str(max(0.5, 6 * row["Percentage"]/100)),
                   arrowsize = str(max(1, 2 * row["Percentage"]/100)))
    #save the plot.
    graph.format = "svg"
    graph.render(name, filepath)

def pathway_wait_in_place(pathway_definitions, pathways_wait_in_place,
                          recurrent_processes):
    #Function to create the list of wait in place processes as an output.
    #----
    #Filter the from and to processess to only those in a wait in place pathway.
    #Get a list of these with no duplicates
    from_col = (pathway_definitions.loc[pathway_definitions['From Process']
                .str.contains('|'.join(pathways_wait_in_place)), 'From Process']
                .drop_duplicates().dropna().to_list())
    to_col = (pathway_definitions.loc[pathway_definitions['To Process']
              .str.contains('|'.join(pathways_wait_in_place)), 'From Process']
              .drop_duplicates().dropna().to_list())
    wip_processes = from_col + to_col
    #Add the list of repeated processes to the wip list
    wip_processes += [process for process in recurrent_processes if
                      any([pathway in process for pathway
                           in pathways_wait_in_place])]
    #Create and return wait in place dataframe
    wip_processes = list(set(wip_processes))
    wait_in_place = pd.DataFrame(
        {'Processes which Wait in Place (Pathway or Recurrent)' : wip_processes,
         'Notes' : [np.nan] * len(wip_processes)})
    return wait_in_place

def get_dfg(events_data, export_event_log_csv,
            export_log_to_csv_after_using_log_converter,
            concept="Event (Pathway)", resource="EventLocation",
            filepath=Path(".")):
    #Function to create the pathway diagram for each stream.
    #tidy data
    events_data = events_data.dropna(subset=["EventTime"])
    events_data = sort_events(events_data)
    events_data.rename(columns={"EventTime": "time:timestamp",
                               "VisitId": "case:concept:name",
                               concept: "concept:name",
                               resource: "org:resource"}, inplace=True)
    #write events log to csv
    if export_event_log_csv:
        events_data.to_csv(str(filepath / "Events Log.csv"), index=False)
    # Convert to log format and save
    log = log_converter.apply(events_data)
    if export_log_to_csv_after_using_log_converter:
        pd.DataFrame(log).to_csv(str(filepath / "Log.csv"), index=False)
    dfg = dfg_discovery.apply(log)
    return dfg

def exclude_patients_with_transitions_below_threshold(events_data,
                                                               threshold):
        #Function to remove patients who have a transition below the specified
        #threshold.
        visit_ids_to_exclude = (events_data
                                .loc[(events_data["EventName"] != "Spawn")
                                    & (events_data["EventName"] != "Removed")
                                    & (events_data["Percentage"] < threshold)
                                    & (events_data["Next Event (Pathway)"]
                                      .notnull()),"VisitId"].unique())
        result = (events_data.loc[~events_data["VisitId"]
                                  .isin(visit_ids_to_exclude)].copy())
        return result

def remove_transitions_below_percentage(pathway_definitions, threshold):
        #Function to remove transitions below a specified threshold.
        filtered_pathway_defs = (pathway_definitions
                                        .loc[pathway_definitions["Percentage"]
                                             > threshold].copy())
        event_pathway_groupby_new_sums = (filtered_pathway_defs
                                          .groupby("From Process")["Percentage"]
                                          .transform("sum"))
        filtered_pathway_defs["Percentage"] = (filtered_pathway_defs
                                               ["Percentage"]
                                            * (100
                                            / event_pathway_groupby_new_sums))
        return filtered_pathway_defs

################################################################################
#----------------------main pathway definitions function-----------------------#
################################################################################

def main_generate_dfg_and_pathway_definitions(directory_path, events_data,
    filepath, include_spawn_end_events, obs_splits, export_event_log_csv,
    export_log_to_csv_after_using_log_converter, pathways_wait_in_place,
    threshold, threshold_exclude, percentage_exclude,
    process_column="Event (Pathway)", split_column=" "):
    #Function to generate and output the pathway definitions file and the
    #DFG/pathway diagrams.
    #----
    #create file directory if it doesn't exist
    filepath.mkdir(exist_ok=True, parents=True)
    #If threshold exclude, filter out events below the threshold.
    if threshold_exclude:
        events_data = exclude_patients_with_transitions_below_threshold(
                      events_data, threshold)
        events_data = add_reset_transitions(events_data)
    #save the events data to csv
    events_data.to_csv(str(filepath / "Events Data.csv"), index=False)
    #get pathway definitions and add in the triage obs events.
    pathway_definitions = generate_pathway_definitions(events_data,
                                                       include_spawn_end_events)
    pathway_definitions = add_obs_repeat_splits(pathway_definitions, obs_splits)
    #Create process recurrence outputs
    process_recurrence_triggers, process_recurrence = create_process_recurrence(
                                                      obs_splits)
    process_recurrence_triggers.to_csv(
        str(filepath / "Process Recurrence Triggers.csv"), index=False)
    process_recurrence.to_csv(
        str(filepath / "Process Recurrence.csv"), index=False)
    #if post processing functions, apply these, then save the pathway
    # definitions to csv and output the transitions plot.
    if percentage_exclude:
        pathway_definitions = remove_transitions_below_percentage(
                              pathway_definitions, threshold)
    pathway_definitions.to_csv(str(filepath / "Pathway Definition.csv"),
                               index=None)
    output_transition_viz(pathway_definitions, filepath, directory_path)
    #Create and save process wait in place.
    wait_in_place = pathway_wait_in_place(pathway_definitions,
                                          pathways_wait_in_place,
                                          [lst[3] for lst in obs_splits])
    wait_in_place.to_csv(filepath/"Process Wait in Place.csv", index=False)
    #Create log files and dfgs
    if split_column is None:
        logs = {}
        logs["All"] = get_dfg(events_data, export_event_log_csv,
                              export_log_to_csv_after_using_log_converter,
                              process_column, filepath=filepath)
    else:
        logs = {}
        for i in events_data[split_column].unique():
            split_filepath = filepath / i
            logs[f"{i}"] = get_dfg(events_data.loc[events_data[split_column]
                                                   == i].copy(),
                                   export_event_log_csv,
                                   export_log_to_csv_after_using_log_converter,
                                   process_column, filepath=split_filepath)
    #save the direct follows graphs.
    for key, dfg in logs.items():
        gviz = dfg_visualizer.apply(dfg)
        dfg_visualizer.save(gviz, filepath / f"{key}.png")
