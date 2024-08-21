import pandas as pd
import os

def print_and_add_str(output_str, added_str):
    print(added_str)
    return output_str + added_str + '\n'

def print_missing(def_str, list, output_str):
    if list:
        output_str = print_and_add_str(output_str, def_str)
        for i in list:
            output_str = print_and_add_str(output_str, i)
    return output_str

output_text = ''
scenario_filepath = 'G:/PerfInfo/Performance Management/PIT Adhocs/2024-2025/HannahP 2425/UEC Adhocs/Frazer Nash UEC/Baseline with resources'
outputs_filepath = 'G:/PerfInfo/Performance Management/PIT Adhocs/2024-2025/HannahP 2425/UEC Adhocs/Frazer Nash UEC/Baseline with resources'
os.chdir(scenario_filepath)

arrival_rates = pd.read_csv('Arrival Rates.csv')
simulation_settings = pd.read_csv('Simulation Settings.csv')
output_text = print_and_add_str(output_text, simulation_settings.to_string())

##################################RESOURCES#####################################
#read in resource csvs
resource_requirement = pd.read_csv('Process Resource Requirement.csv').add_suffix(' Res Req')
resource_rota = pd.read_csv('Resource Rota.csv').add_suffix(' Res Rota')

#Split up resource requirement column to break up the different staff members
#being used for each task.  Explode this so one row per staff member, and split
#up the staff memeber name and the number of that staff required.
resource_requirement['Requirement Res Req'] = (resource_requirement['Requirement Res Req']
                                       .str.split(', '))
resource_requirement = resource_requirement.explode(column='Requirement Res Req')
resource_requirement[['Number of Resource Res Req',
                      'Requirement Res Req']] = (resource_requirement['Requirement Res Req']
                                         .str.split(r'(?<=\d)\s', regex=True,
                                                    expand=True))

#merge this onto the resource rota
resources = (resource_requirement.merge(resource_rota,
                                        left_on='Requirement Res Req', 
                                        right_on='Resource Type Res Rota',
                                        how='outer'))

output_text = print_and_add_str(output_text, '-----------------RESOURCES-----------------')
#print the staff that aren't in use
unused_staff = resources.loc[resources['Requirement Res Req'].isna(),
                        'Resource Type Res Rota'].drop_duplicates().dropna().to_list()
output_text = print_missing('--Unused staff in Resource Rota are:', unused_staff, output_text)

#print the staff that aren't scheduled
unscheduled_staff = resources.loc[resources['Resource Type Res Rota'].isna(),
                            'Requirement Res Req'].drop_duplicates().dropna().to_list()
output_text = print_missing("--Staff required to complete processess that aren't in the rota:",
              unscheduled_staff, output_text)

##################################LOCATIONS#####################################
#read in location csvs
process_locations = pd.read_csv('Process Locations.csv').add_suffix(' Proc Loc')
location_capacities = pd.read_csv('Location Capacities.csv').add_suffix(' Loc Cap')
location_opening_hours = pd.read_csv('Location Opening Hours.csv').add_suffix(' Loc Open')

#merge into one locations df
locations = (process_locations.merge(location_capacities,
                                     left_on='Location Proc Loc',
                                     right_on='Location Loc Cap',
                                     how='outer')
                              .merge(location_opening_hours,
                                       left_on='Location Loc Cap',
                                       right_on='Location Loc Open',
                                       how='outer'))

output_text = print_and_add_str(output_text, '-----------------LOCATIONS-----------------')
#Print the process locations with capacities missing
proc_loc_no_cap = (locations.loc[locations['Location Loc Cap'].isna(),
                                 'Location Proc Loc'].drop_duplicates()
                                 .dropna().to_list())
output_text = print_missing('--Process locations with missing location capacities:',
              proc_loc_no_cap, output_text)
#Print location capacities that don't match process locations
cap_no_proc_loc = (locations.loc[locations['Location Proc Loc'].isna(),
                                 'Location Loc Cap'].drop_duplicates()
                                 .dropna().to_list())
output_text = print_missing("--Locations with capacities that don't match to a process location:",
              cap_no_proc_loc, output_text)

#Print the locations that don't have opening hours
no_open_hours = locations.loc[locations['Location Loc Open'].isna(),
                    'Location Proc Loc'].drop_duplicates().dropna().to_list()
output_text = print_missing("--Locations that don't have opening hours set:", no_open_hours, output_text)

#print the locations that have opening hours but no capacity
no_capacity = locations.loc[locations['Location Proc Loc'].isna(),
                    'Location Loc Open'].drop_duplicates().dropna().to_list()
output_text = print_missing("--Locations that don't have capacity data:", no_capacity, output_text)



#############################RECURRENT PROCESSES################################
#read in csvs
process_recurrence_triggers = pd.read_csv('Process Recurrence Triggers.csv').add_suffix(' Rec Trig')
process_recurrence = pd.read_csv('Process Recurrence.csv').add_suffix(' Pro Rec')
#merge into one df
recurrence = (process_recurrence_triggers.merge(process_recurrence,
              left_on='Recurrent Process (Not In Pathway) Rec Trig', 
              right_on='Recurrent Process Pro Rec', how='outer'))

output_text = print_and_add_str(output_text, '---------------RECURRENT PROCESSES---------------')
output_text = print_and_add_str(output_text, 'NOTE: This only checks matches between the two recurrent process files.')
output_text = print_and_add_str(output_text, '      Comparison to pathway file done later.')
#Print the process in pathway that doesn't match to a recurrent process
no_rec_proc = (recurrence.loc[recurrence['Recurrent Process Pro Rec'].isna(),
                              'Recurrent Process (Not In Pathway) Rec Trig']
                              .drop_duplicates().dropna().to_list())
output_text = print_missing('--Process in pathway trigger without matching recurrent processes:',
              no_rec_proc, output_text)

#Print recurrent processes that don't match to the pathway process
no_path_proc = (recurrence.loc[
                recurrence['Recurrent Process (Not In Pathway) Rec Trig'].isna(),
                'Recurrent Process Pro Rec'].drop_duplicates().dropna().to_list())
output_text = print_missing('--Process in processes recurrence but not in trigger file:',
              no_path_proc, output_text)

##################################PROCESSES#####################################
#read in csvs
process_durations = pd.read_csv('Process Durations.csv').add_suffix(' Proc Dur')
process_wait_in_place = pd.read_csv('Process Wait in Place.csv').add_suffix(' Proc Wait')
arrival_rates = pd.read_csv('Arrival Rates.csv').add_suffix(' Arr Rat')

#merge into one df
processes = (process_durations.merge(process_wait_in_place,
                                      left_on='Process (Pathway and Recurrent) Proc Dur',
                                      right_on='Processes which Wait in Place (Pathway or Recurrent) Proc Wait',
                                      how='outer')
                              .merge(arrival_rates, left_on='Process (Pathway and Recurrent) Proc Dur',
                                     right_on='InitialProcess Arr Rat', how='outer'))

#Print the processess in the durations file that aren't in the locations file
output_text = print_and_add_str(output_text, '---------------PROCESSES DURATIONS AND WIP---------------')
#Print the processess in the durations file that don't match wait in place
no_wip = (processes
          .loc[processes['Processes which Wait in Place (Pathway or Recurrent) Proc Wait'].isna(),
               'Process (Pathway and Recurrent) Proc Dur']
               .drop_duplicates().dropna().to_list())
output_text = print_missing("--Process Durations that aren't Wait in Place (No Issue):", no_wip, output_text)

#Print the processes in the Wait in Place file which don't match the durations
#file
wip_no_dur = (processes.loc[processes['Process (Pathway and Recurrent) Proc Dur'].isna(),
                            'Processes which Wait in Place (Pathway or Recurrent) Proc Wait']
                            .drop_duplicates().dropna().to_list())
output_text = print_missing("--Wait in Place processes that don't match the Process Durations:",
              wip_no_dur, output_text)

#Print the arrival processes that don't map to a duration process
arr_proc_no_dur = (processes.loc[(processes['Process (Pathway and Recurrent) Proc Dur'].isna())
                                 & ~(processes['InitialProcess Arr Rat'].isna()),
                                 'InitialProcess Arr Rat']
                                 .drop_duplicates().dropna().to_list())
output_text = print_missing("--Arrival processes that don't match the Process Durations:",
              arr_proc_no_dur, output_text)



##############################PATHWAY DEFINITION################################
#read in csv
output_text = print_and_add_str(output_text, '------------------------------------------------------')
output_text = print_and_add_str(output_text, '------------MERGING ON PATHWAY DEFINITIONS------------')
output_text = print_and_add_str(output_text, '------------------------------------------------------')
output_text = print_and_add_str(output_text, ' ')
pathway_definition = pd.read_csv('Pathway Definition.csv').add_suffix(' Path Def')

def pathway_def_merge(pathway_def, from_or_to, output_text):
    ####Merge onto other dataframes
    output_text = print_and_add_str(output_text, f'------------{from_or_to} PROCESSES------------')
    output_text = print_and_add_str(output_text, ' ')
    pathway_proc_col = f'{from_or_to} Process Path Def'
    pathway_def = (pathway_def.merge(processes, left_on=pathway_proc_col,
                                  right_on='Process (Pathway and Recurrent) Proc Dur',
                                  how='outer')
                              .merge(recurrence, left_on=pathway_proc_col,
                                     right_on='Trigger Process (In Pathway) Rec Trig',
                                     how='outer')
                              .merge(resources, left_on=pathway_proc_col,
                                     right_on='Process (Pathway or Recurrent) Res Req',
                                     how='outer')
                              .merge(locations, left_on=pathway_proc_col,
                                     right_on='Process Proc Loc', how='outer'))

    pathway_def_simp = pathway_def[[pathway_proc_col,
                                  'Process (Pathway and Recurrent) Proc Dur',
                                  'Trigger Process (In Pathway) Rec Trig',
                                  'Recurrent Process (Not In Pathway) Rec Trig',
                                  'Recurrent Process Pro Rec',
                                  'InitialProcess Arr Rat',
                                  'Process (Pathway or Recurrent) Res Req',
                                  'Requirement Res Req',
                                  'Resource Type Res Rota', 'Process Proc Loc',
                                  'Location Proc Loc', 'Location Loc Cap',
                                  'Location Loc Open']]

    #get rows where pathway column is nan
    pathway_nas = pathway_def[pathway_proc_col].isna()

    output_text = print_and_add_str(output_text, '---------------PROCESSES DURATIONS---------------')
    #get rows where process durations is nan
    proc_dur_nas = pathway_def['Process (Pathway and Recurrent) Proc Dur'].isna()
    #print processess in the pathway that are missing durations
    proc_no_dur = (pathway_def.loc[(proc_dur_nas) & ~(pathway_nas),
                                   pathway_proc_col]
                                   .drop_duplicates().dropna().to_list())
    output_text = print_missing('--Pathway Processes with missing durations:', proc_no_dur, output_text)
    #Print durations that aren't in pathway processes
    dur_no_proc = (pathway_def.loc[(pathway_nas) & ~(proc_dur_nas), 
                                   'Process (Pathway and Recurrent) Proc Dur']
                                   .drop_duplicates().dropna().to_list())
    output_text = print_missing('--Process durations not in pathway definition:', dur_no_proc, output_text)


    output_text = print_and_add_str(output_text, '---------------RECURRENT PROCESSES---------------')
    #get rows where recurrent triggers are nan
    rec_trig_nas = pathway_def['Trigger Process (In Pathway) Rec Trig'].isna()
    #print recurrent processes that aren't in the pathway definitions.
    missing_rec = (pathway_def.loc[(pathway_nas) & ~(rec_trig_nas),
                                   'Trigger Process (In Pathway) Rec Trig']
                                   .drop_duplicates().dropna().to_list())
    output_text = print_missing('--Recurrent Processes missing from pathway:', missing_rec, output_text)

    output_text = print_and_add_str(output_text, '---------------RESOURCES---------------')
    #get rows where resource requirement is na
    res_req_nas = pathway_def['Process (Pathway or Recurrent) Res Req'].isna()
    #Print processes in the pathway that don't have resources
    path_no_res = (pathway_def.loc[(res_req_nas) & ~(pathway_nas),
                                   pathway_proc_col]
                                   .drop_duplicates().dropna().to_list())
    output_text = print_missing("--Pathway processes with no staff requirements:", path_no_res, output_text)
    #Print resource processes which aren't in the pathway definition
    res_no_path = (pathway_def.loc[(pathway_nas) & ~(res_req_nas),
                                   'Process (Pathway or Recurrent) Res Req']
                                   .drop_duplicates().dropna().to_list())
    output_text = print_missing("--Resource processes that are missing from the pathway definition:",
                  res_no_path, output_text)
    
    output_text = print_and_add_str(output_text, '---------------LOCATIONS---------------')
    #get rows where Location is na
    proc_loc_nas = pathway_def['Process Proc Loc'].isna()

    #Get the pathway definitions that are missing process locations
    path_no_proc_loc = (pathway_def.loc[(proc_loc_nas) & ~(pathway_nas),
                                        pathway_proc_col]
                                        .drop_duplicates().dropna().to_list())
    output_text = print_missing("--Pathway processes that don't match with a process location:",
                  path_no_proc_loc, output_text)
    #print the processes locations that aren't in pathway definition
    proc_loc_no_path = (pathway_def.loc[(pathway_nas) & ~(proc_loc_nas),
                                        'Process Proc Loc']
                                        .drop_duplicates().dropna().to_list())
    output_text = print_missing('--Process Locations that are missing from the pathway definition:',
                  proc_loc_no_path, output_text)
    
    return pathway_def, pathway_def_simp, output_text

from_pathway_def, from_pathway_simp, output_text = pathway_def_merge(pathway_definition.copy(), 'From', output_text)
to_pathway_def, to_pathway_simp, output_text = pathway_def_merge(pathway_definition.copy(), 'To', output_text)

#Change to output filepath (create if doens't exist)
validation_filepath = outputs_filepath + '/validation'
if not os.path.exists(validation_filepath):
    os.mkdir(validation_filepath)
os.chdir(validation_filepath)

#Export smaller combined data frames
processes.to_csv('processes validaiton.csv', index=False)
locations.to_csv('locations validation.csv', index=False)
resources.to_csv('resources validation.csv', index=False)
recurrence.to_csv('recurrence validation.csv', index=False)
#export from and to combined data frames
from_pathway_def.to_csv('combined from pathway validation.csv', index=False)
from_pathway_simp.to_csv('combined from pathway validaton simple.csv', index=False)
to_pathway_def.to_csv('combined to pathway validation.csv', index=False)
to_pathway_simp.to_csv('combined to pathway validation simple.csv', index=False)

with open('Output.txt', 'w', encoding='utf-8') as f:
    f.write(output_text)