import pandas as pd
import numpy as np
from datetime import datetime

############################FILEPATHS############################
output_path = "G:/PerfInfo/Performance Management/PIT Adhocs/2024-2025/HannahP 2425/UEC Adhocs/Frazer Nash UEC/Code Scenarios/" + datetime.today().strftime('%Y-%m-%d')
other_input_filepath = "G:/PerfInfo/Performance Management/PIT Adhocs/2024-2025/HannahP 2425/UEC Adhocs/Frazer Nash UEC/Baseline with resources"

#######################THRESHOLD VARIABLES#######################
collapse_diagnostics_within_time = pd.Timedelta("5m")
transition_threshold = 2
quantile_threshold = 1.0
repeat_time_threshold = 10

#############################BOOLS###############################
remove_duplicate_staffid = False
remove_duplicate_location = False
include_obs_data = True
include_diag_data = True
include_admission_data = True
export_event_log_csv = False
export_log_to_csv_after_using_log_converter = False
keep_last_location = True
include_spawn_end_events = False
plots = True
print_validation = False

###########################LISTS/DICTS###########################
pathways =  ["Minors", "Ambulatory", "Majors", "Resus"]
additional_filenames = ["Arrival Rates", "Location Opening Hours",
                        "Simulation Settings", "Location Capacities",
                        "Resource Rota", "Process Resource Requirement",
                        "Process Locations"]
event_names_to_exclude_for_repetition = ["Triaged", "Discharged", "Booked In",
                                         "Ambulance Arrival", "Walk-In",
                                         "Admitted - Other Derriford Ward",
                                         "Admitted - MAU", "Admitted - SDEC"]
#where_duration_should_be_0 = ["Walk-In",
 #                             "Admitted - Other Derriford Ward",
  #                            "Admitted - MAU", "Admitted - SDEC", "Discharged",
   #                           "Admitted",
    #                          "Triaged - Kickoff 60 min Obs",
     #                         "Triaged - Kickoff 30 min Obs", "Spawn", "Removed"]
imaging_events = ["Radiology", "CT", "MRI", "Ultrasound"]
locations_pathway_map = {"Ambulance": "Majors",
                         "Ambulatory Cubicles": "Ambulatory",
                         "Ambulatory Waiting Area": "Ambulatory",
                         "Majors Corridor": "Majors",
                         "Majors Cubicles": "Majors",
                         "Minors": "Minors",
                         "Resus": "Resus"}
admitted_map = {'Admitted - MAU':'Admitted',
                'Admitted - Other Derriford Ward':'Admitted',
                'Admitted - SDEC':'Admitted'}
pathways_wait_in_place = ["Majors", "Resus"]
excluded_event_names = ["Clinically Ready to Proceed"]
locations_to_drop = ["Paediatrics", "Plym"]
location_opening_hours = {'Location':[], 'Day of Week':[], 'Start Time':[],
                          'End Time':[], 'Notes':[]}
#Process Durations that are manually added
# process : [mean, std, min, max]
add_process_durs = {#'Admitted':[0, 0, 0, 0],
                    #'Discharged':[0, 0, 0, 0],
                    'Triaged':[10, 6, 1, 40],
                    'Booked In':[3, 2, 1, 10],
                    'CT':[20, 10, 10, 40],
                    'Radiology':[10, 5, 5, 30],
                    'Obs 60 min (Ambulatory)':[10, 1, 5, 20],
                    'Obs 60 min (Resus)':[10, 1, 5, 20],
                    'Obs 60 min (Majors)':[10, 1, 5, 20],
                    'Obs 30 min (Resus)':[10, 1, 5, 20],
                    'Obs 30 min (Majors)':[10, 1, 5, 20],
                    'Misc Assessment' : [10, 5, 5, 30]}
proc_durs_0 = ['Admitted (Minors)',
               'Admitted (Ambulatory)',
               'Admitted (Majors)',
               'Admitted (Resus)',
               'Discharged (Minors)',
               'Discharged (Ambulatory)',
               'Discharged (Majors)',
               'Discharged (Resus)',
               'Walk-In (Ambulatory)',
               'Walk-In (Majors)',
               'Walk-In (Minors)',
               'Walk-In (Resus)',
               'Triaged - Kickoff 60 min Obs (Majors)',
               'Triaged - Kickoff 30 min Obs (Majors)',
               'Triaged - Kickoff 60 min Obs (Resus)',
               'Triaged - Kickoff 30 min Obs (Resus)']
#Put 0 time processes into dictionary to be added
for process in proc_durs_0:
    add_process_durs[process]  = [0, np.nan, np.nan, np.nan]

manual_process_timings = set([proc.split(' (')[0] for proc in add_process_durs.keys()])
#list of the new events to add after triage to kick off repeated obs and their
#probabilities.
#[(From Event, To Event, Probability, Recurrent Process)]
#From Event is nan if no kickoff event required to start the repeated process.
obs_splits = [(np.nan, 'Triaged (Ambulatory)', 100, 'Obs 60 min (Ambulatory)'),
              ('Triaged (Majors)', 'Triaged - Kickoff 60 min Obs (Majors)', 85,
               'Obs 60 min (Majors)'),
              ('Triaged (Majors)', 'Triaged - Kickoff 30 min Obs (Majors)', 15,
               'Obs 30 min (Majors)'),
              ('Triaged (Resus)', 'Triaged - Kickoff 60 min Obs (Resus)', 60,
               'Obs 60 min (Resus)'),
              ('Triaged (Resus)', 'Triaged - Kickoff 30 min Obs (Resus)', 40,
               'Obs 30 min (Resus)')]

########################PROCESS ORDER########################
natural_order_for_processes = { "Spawn": 0,
                                "Walk-In": 1,
                                "Ambulance Arrival": 1,
                                "Booked In": 2,
                                "Triaged": 3,
                                "Nurse Assessment": 4,
                                "Misc Assessment" : 5,
                                "Bloods/ECG" : 6,
                                "Treatment" : 7,
                                "Seen By Clinician/Treated": 8,
                                "Clerked": 9,
                                "Observations": 10,
                                #"Imaging": 8,
                                "Radiology" : 11,
                                "CT" : 12,
                                "MRI" : 13,
                                "Ultrasound" : 14,
                                "Laboratory": 15,
                                "Specialty Reviewed": 16,
                                "Senior Reviewed": 17,
                                "Clinically Ready to Proceed": 18,
                                "Decision to Admit": 19,
                                "Discharged": 20,
                                "Wait for Bed - Admitted - MAU": 20,
                                "Wait for Bed - Admitted - Other Derriford Ward": 20,
                                "Wait for Bed - Admitted - SDEC": 20,
                                "Admitted": 21,
                                "Removed": 22}
