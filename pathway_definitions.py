from utils import sort_events
from pathlib import Path
from pm4py.objects.conversion.log import converter as log_converter  # type: ignore
from pm4py.algo.discovery.dfg import algorithm as dfg_discovery  # type: ignore
from pm4py.visualization.dfg import visualizer as dfg_visualizer  # type: ignore
from config import SPAWN, REMOVED
from graphviz import Digraph
import pandas as pd


def add_reset_transitions(events_data):
    """
    Args:
        events_data (pd.DataFrame): clensed events dataframe.

    Returns:
        pd.DataFrame: clensed events dataframe with the count and perctage of
        each transition.
    """
    transitions = events_data.copy()
    transitions = transitions.drop(["Count", "Total", "Percentage"], axis=1,
                                   errors="ignore", inplace=False)

    # Calculate Next Event for each Patient
    transitions["Next Event (Pathway)"] = (transitions.groupby("VisitId")
                                           ["Event (Pathway)"].shift(-1))

    # Calculate the count of patients that make each transition
    count_event_pairs = (transitions.groupby(["Event (Pathway)",
                                              "Next Event (Pathway)"],
                                              as_index=False)
                                              .agg({"VisitId": "count"})
                                              .rename(columns=
                                                      {"VisitId": "Count"})
                                              .drop_duplicates())
    #calculate the total number of patients that have each event
    count_events = (transitions.loc[~transitions["Next Event (Pathway)"]
                                    .isnull()].groupby(["Event (Pathway)"],
                                                       as_index=False)
                                    .agg({"VisitId": "count"})
                                    .rename(columns={"VisitId": "Total"})
                                    .drop_duplicates())
    #merge counts onto events data
    transitions = (transitions.merge(count_event_pairs,
                                     on=["Event (Pathway)",
                                         "Next Event (Pathway)"],
                                         how="left")
                              .merge(count_events, on=["Event (Pathway)"],
                                     how="left"))
    #calculate the percentage of patients that make each transition from each
    #event
    transitions["Percentage"] = 100 * transitions["Count"] / transitions["Total"]

    return transitions


def remove_transitions_below_percentage_in_pathway_definitions(threshold):
    """
    Args:
        threshold (int): percentage threshold to remove transitions.
    """
    def remove_lines_from_pathway_definition(pathway_definitions):
        """
        Args:
            pathway_definitions (pd.DataFrame): Dataframe of pathway definitions.

        Returns:
            pd.DataFrame: Dataframe with transitions below threshold removed.
        """
        filtered_pathway_definitions = (pathway_definitions
                                        .loc[pathway_definitions["Percentage"]
                                             > threshold].copy())
        event_pathway_groupby_new_sums = (filtered_pathway_definitions
                                          .groupby("From Process")["Percentage"]
                                          .transform("sum"))
        filtered_pathway_definitions["Percentage"] = (filtered_pathway_definitions["Percentage"]
                                                      * (100 / event_pathway_groupby_new_sums))
        return filtered_pathway_definitions

    return remove_lines_from_pathway_definition


def generate_and_output_pathway_definitions(events_data,
                                            include_spawn_end_events):
    """
    Args:
        events_data (pd.DataFrame): clesed events dataframe.
        include_spawn_end_events (bool): flage to include spawn and end events.

    Returns:
        pd.DataFrame: dataframe of the pathway definitions.
    """
    #If spawn and end events were included, remove them here for next steps.
    if include_spawn_end_events:
        events_data = events_data.loc[~events_data["EventName"]
                                      .isin([SPAWN, REMOVED])].copy()

    #Remove duplicates from to events, remove unecessary columns and rename.
    pathway_definition = (events_data.drop_duplicates(
                         ["Event (Pathway)", "Next Event (Pathway)"])
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
    """
    Args:
        pathway_definition (pd.DataFrame): Pathway definition dataframe
        obs_splits (list(tup(str, str, int)): list of the new events to add
        after triage to kick off repeated obs and their probabilities.

    Returns:
        pd.DataFrame: The pathway definition dataframe with the new events to add
        after triage to kick off repeated obs added
    """
    #Filter the pathway to the triage events that need the extra event to
    #trigger different obs repetitions. Create df of these and remove them from
    #the pathway definitions.
    triage_events_for_obs = set([i[0] for i in obs_splits])
    triage_events_mask = pathway_definition['From Process'].isin(triage_events_for_obs)
    triage_events = pathway_definition.loc[triage_events_mask].copy()
    pathway_definition = pathway_definition.loc[~triage_events_mask].copy()

    #Loop through the list of splits, and create a new dataframe with triaged
    #events going to the triaged obs events with their probabilities.
    triage_obs = []
    for triage_values in obs_splits:
        triage_obs.append(
            [triage_values[0], triage_values[1], None, float(triage_values[2]),
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
    pathway_definition = pd.concat([triage_obs_events, new_links, pathway_definition])

    return pathway_definition


def get_dfg(event_data, export_event_log_csv,
            export_log_to_csv_after_using_log_converter,
            concept="Event (Pathway)", resource="EventLocation",
            filepath=Path(".")):
    """
    Args:
        event_data (pd.DataFrame): clensed events dataframe.
        export_event_log_csv (bool): flag to import event log csv.
        export_log_to_csv_after_using_log_converter (bool): ?????????????
        concept (str, optional): ????????. Defaults to "Event (Pathway)".
        resource (str, optional): string of the column to use.
        Defaults to "EventLocation".
        filepath (Path, optional):filepath to save results. Defaults to Path(".").

    Returns:
        dfg: Directly Follows Graph is returned.
    """
    #tidy data
    event_data = event_data.dropna(subset=["EventTime"])
    event_data = sort_events(event_data)
    event_data.rename(columns={"EventTime": "time:timestamp",
                               "VisitId": "case:concept:name",
                               concept: "concept:name",
                               resource: "org:resource"},inplace=True)
    #write events log to csv
    if export_event_log_csv:
        event_data.to_csv(str(filepath / "Events Log.csv"), index=False)

    # Convert to log format and save
    log = log_converter.apply(event_data)
    if export_log_to_csv_after_using_log_converter:
        pd.DataFrame(log).to_csv(str(filepath / "Log.csv"), index=False)

    dfg = dfg_discovery.apply(log)
    return dfg


def output_transition_viz(pathway_definition, filepath, name):
    """
    Args:
        pathway_definition (pd.DataFrame): pathway definitions dataframe.
        filepath (str): filepath to save the visualisation to.
        name (str): name of the folder being saved to.
    """

    def get_colour(node):
        """
        Args:
            node (str): the nodes of the plot.

        Returns:
            str: the output colour of each node.
        """
        if 'Arrival' in node:
            return 'lightcoral'
        if 'Walk-In' in node:
            return 'lightcoral'
        
        if 'Admitted' in node and not 'Wait' in node:
            return 'lightgreen'
        if 'Discharged' in node and not 'Wait' in node:
            return 'lightgreen'
        return 'aliceblue'
    
    def get_shape(node : str) -> str:
        """
        Args:
            node (str): the nodes of the plot.

        Returns:
            str: the output shape of each node.
        """
        if 'Arrival' in node:
            return 'invhouse'
        if 'Walk-In' in node:
            return 'invhouse'
        
        if 'Admitted' in node and not 'Wait' in node:
            return 'house'
        if 'Discharged' in node and not 'Wait' in node:
            return 'house'
        return 'oval'

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


def generate_and_output_dfg_and_pathway_definition(directory_path, events_data,
    filepath, include_spawn_end_events, obs_splits, export_event_log_csv,
    export_log_to_csv_after_using_log_converter, location_capacity_data,
    process_location_data, event_names_based_process_requirements,
    pathways_wait_in_place, process_column="Event (Pathway)",
    split_column=" ", filterFuncs=None,
    post_processing_functions_of_pathway_definitions=None):
    """
    Args:
        directory_path (str): name of folder to create or populate.
        events_data (pd.DataFrame): clensed events dataframe.
        filepath(str): full filepath to folder to populate.
        include_spawn_end_events (bool): flag to include spawn events.
        obs_splits(list(tuple(str, str, int))): list of the new events to add
        after triage to kick off repeated obs and their probabilities.
        export_event_log_csv (bool): flag to export event log csv.
        export_log_to_csv_after_using_log_converter (bool): flag to export
        log csv.
        location_capacity_data (dict [str, str | int]): dictionary of the
        location capacities.
        process_location_data (list[tuple[str, str, int]]): list of tuples of
        the process locations.
        event_names_based_process_requirements (list[tuple[str, str, int]]): 
        list of tuples of the resources required for each process.
        pathways_wait_in_place (list[str]): list of pathways where patient
        waits in place.
        process_column (str): the column name of the events,
        default="Event (Pathway)".
        split_column (str): column to split logs by if required. default=" "
        filterFuncs (): list of filter functions if required default=None
        post_processing_functions_of_pathway_definitions (): list of functions
        for postprocessing of process definitions if required. default=None
    """
    #create file directory if it doesn't exist
    filepath.mkdir(exist_ok=True, parents=True)

    if filterFuncs is not None:
        for filter in filterFuncs:
            events_data = filter(events_data)
            events_data = add_reset_transitions(events_data)

    #save the events data to csv
    events_data.to_csv(str(filepath / "Events Data.csv"), index=False)

    #get pathway definitions and add in the triage obs events.
    pathway_definitions_unfiltered = generate_and_output_pathway_definitions(
                                     events_data, include_spawn_end_events)
    pathway_definitions_unfiltered = add_obs_repeat_splits(pathway_definitions_unfiltered,
                                                           obs_splits)

    #if post processing functions, apply these, then save the pathway
    # definitions to csv and output the transitions plot.
    if post_processing_functions_of_pathway_definitions is not None:

        filtered_pathway_definitions = post_processing_functions_of_pathway_definitions(
                                       pathway_definitions_unfiltered)
        filtered_pathway_definitions.to_csv(filepath / "Pathway Definition.csv",
                                            index=None)
        output_transition_viz(filtered_pathway_definitions, filepath,
                              directory_path)

    else:
        pathway_definitions_unfiltered.to_csv(filepath/"Pathway Definition.csv",
                                              index=False)
        output_transition_viz(pathway_definitions_unfiltered, filepath,
                              directory_path)

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
 #           split_filepath.mkdir(exist_ok=True, parents=True)
            logs[f"{i}"] = get_dfg(events_data.loc[events_data[split_column]
                                                   == i].copy(),
                                   export_event_log_csv,
                                   export_log_to_csv_after_using_log_converter,
                                   process_column, filepath=split_filepath)
    #save the direct follows graphs.
    for key, dfg in logs.items():
        gviz = dfg_visualizer.apply(dfg)
        dfg_visualizer.save(gviz, filepath / f"{key}.png")

