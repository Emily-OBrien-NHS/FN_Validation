import pandas as pd
import numpy as np

#######################THRESHOLD VARIABLES#######################
collapse_diagnostics_within_time = pd.Timedelta("5m")
quantile_threshold = 0.97
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

###########################LISTS/DICTS###########################
event_names_to_exclude_for_repetition = ["Triaged", "Discharged", "Booked In",
                                         "Ambulance Arrival", "Walk-In",
                                         "Admitted - Other Derriford Ward",
                                         "Admitted - MAU", "Admitted - SDEC"]
where_duration_should_be_0 = ["Wait for Bed", "Walk-In",
                              "Admitted - Other Derriford Ward",
                              "Admitted - MAU", "Admitted - SDEC", "Discharged"]
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
