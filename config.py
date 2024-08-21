"""
This module defines configuration parameters for analysis purposes.

The module includes constants, time variables, boolean flags, excluded objects,
lists for additional outputs, and mappings used throughout the analysis process.

Configuration parameters:
- Constants: Predefined string constants for event names and pathways.
- Time Variables: Time durations used in analysis.
- Boolean Flags: Toggles for including various types of data in the analysis.
- Lists: Collections of strings representing excluded event names, locations, and pathways.
- Dictionaries: Mappings for natural order of processes, location pathways, and capacity data.

"""
import pandas as pd

#######################THRESHOLD VARIABLES#######################
duration_processes_quantile_threshold = 1
MAX_DIFF_MINUTES_FOR_DURATIONS = 840
collapse_diagnostics_rows_within_time_of = pd.Timedelta("5m")
repeat_time_threshold = 10

#######################STRINGS#######################
#Nodes
SPAWN = "Spawn"
REMOVED = "Removed"
WALK_IN = "Walk-In"
ALL_PATHWAYS = "All"
WAITING_FOR_BED = "Wait for bed"
#Staff
RECEPTIONIST = "Receptionist"
CONSULTANT = "Consultant"
MIDDLE = "Middle Grade"
JUNIOR = "Junior Doctor"
NURSE = "Nurse"
ENP = "ENP"
ACP = "ACP"
#Other
PROCESS_REQUIREMENTS_COLUMN = "Event (Pathway)"

#######################BOOLS#######################
remove_duplicate_staffid = False
remove_duplicate_location = False
include_obs_data = False
include_diag_data = False
include_admission_data = True
export_event_log_csv = False
export_log_to_csv_after_using_log_converter = False
keep_last_location = True
include_spawn_end_events = False
plots = False

#######################LISTS/DICTS#######################
event_names_to_exclude_for_repetition = ["Triaged", "Discharged", "Booked In",
                                         "Ambulance Arrival", WALK_IN,
                                         "Admitted - Other Derriford Ward",
                                         "Admitted - MAU", "Admitted - SDEC"]
where_duration_should_be_0 = [WAITING_FOR_BED, WALK_IN,
                              "Admitted - Other Derriford Ward",
                              "Admitted - MAU", "Admitted - SDEC", "Discharged"]
locations_pathway_map = {"Ambulance": "Majors",
                         "Ambulatory Cubicles": "Ambulatory",
                         "Ambulatory Waiting Area": "Ambulatory",
                         "Majors Corridor": "Majors",
                         "Majors Cubicles": "Majors",
                         "Minors": "Minors",
                         "Resus": "Resus"}
pathways_wait_in_place = ["Majors", "Resus"]
excluded_event_names = ["Nursing Assessment", "Clinically Ready to Proceed"]
locations_to_drop = ["Paediatrics", "Plym"]
pathways = ["Majors", "Ambulatory", "Minors", "Resus"]
#list of the new events to add after triage to kick off repeated obs and their
#probabilities.
obs_splits = [('Triaged (Majors)', 'Triaged 60 min Obs (Majors)', 85),
                ('Triaged (Majors)', 'Triaged 30 min Obs (Majors)', 15),
                ('Triaged (Resus)', 'Triaged 60 min Obs (Resus)', 60),
                ('Triaged (Resus)', 'Triaged 30 min Obs (Resus)', 40)]

#######################PROCESS ORDER#######################
natural_order_for_processes = { SPAWN: 0,
                                WALK_IN: 1,
                                "Ambulance Arrival": 1,
                                "Booked In": 2,
                                "Triaged": 3,
                                "Nursing Assessment": 4,
                                "Seen By Clinician/Treated": 5,
                                "Clerked": 6,
                                "Observations": 7,
                                "Imaging": 8,
                                "Laboratory": 9,
                                "Specialty Reviewed": 10,
                                "Senior Reviewed": 11,
                                "Clinically Ready to Proceed": 12,
                                "Decision to Admit": 13,
                                "Discharged": 14,
                                WAITING_FOR_BED + " - (Admitted - MAU)": 14,
                                WAITING_FOR_BED + " - (Admitted - Other Derriford Ward)": 14,
                                WAITING_FOR_BED + " - (Admitted - SDEC)": 14,
                                "Admitted - Other Derriford Ward": 15,
                                "Admitted - MAU": 15,
                                "Admitted - SDEC": 15,
                                REMOVED: 16}

####################PROCESS RESOURCE REQUIREMENTS DATA####################
event_names_based_process_requirements_event_staff_mapping_and_priority = [
    ("Clerked", f"1 {CONSULTANT}", 4),
    ("Clerked", f"1 {MIDDLE}", 3),
    ("Clerked", f"1 {JUNIOR}", 2),
    ("Clerked", f"1 {ACP}", 1),
    ("Triaged", f"1 {NURSE}", 1),
    ("Triaged", f"1 {ENP}", 2),
    ("Triaged", f"1 {ACP}", 3),
    ("Senior Reviewed", f"1 {CONSULTANT}", 2),
    ("Senior Reviewed", f"1 {MIDDLE}", 1),
    ("Decision to Admit", f"1 {CONSULTANT}", 3),
    ("Decision to Admit", f"1 {MIDDLE}", 2),
    ("Decision to Admit", f"1 {ACP}", 1),
    ("Seen By Clinician/Treated", f"1 {CONSULTANT}", 6),
    ("Seen By Clinician/Treated", f"1 {MIDDLE}", 5),
    ("Seen By Clinician/Treated", f"1 {JUNIOR}", 4),
    ("Seen By Clinician/Treated", f"1 {ACP}", 3),
    ("Seen By Clinician/Treated", f"1 {ENP}", 2),
    ("Seen By Clinician/Treated", f"1 {NURSE}", 1),
    ("Discharged", f"1 {CONSULTANT}", 6),
    ("Discharged", f"1 {MIDDLE}", 5),
    ("Discharged", f"1 {JUNIOR}", 4),
    ("Discharged", f"1 {ACP}", 3),
    ("Discharged", f"1 {ENP}", 2),
    ("Discharged", f"1 {NURSE}", 1),
    ("Admitted", f"1 {CONSULTANT}", 6),
    ("Admitted", f"1 {MIDDLE}", 5),
    ("Admitted", f"1 {JUNIOR}", 4),
    ("Admitted", f"1 {ACP}", 3),
    ("Admitted", f"1 {ENP}", 2),
    ("Admitted", f"1 {NURSE}", 1)]

event_names_based_process_requirements = []
for i in event_names_based_process_requirements_event_staff_mapping_and_priority:
    for pathway in pathways:
        event_names_based_process_requirements.append(
            (f"{i[0]} ({pathway})", f"{i[1]}", i[2]))

#######################PROCESS LOCATIONS DATA#######################
process_location_data = [
    ("Ambulance Arrival (Majors)", "Ambulance", 1),
    ("Ambulance Arrival (Resus)", "Ambulance", 1),
    ("Ambulance Arrival (Ambulatory)", "Ambulance", 1),
    ("Ambulance Arrival (Minors)", "Ambulance", 1),
    ("Booked In (Ambulatory)", "Reception", 1),
    ("Booked In (Minors)", "Reception", 1),
    ("Triaged (Resus)", "RAT Room", 1),
    ("Triaged (Resus)", "Resus Bays", 2),
    ("Triaged (Resus)", "Ambulance", 3),
    ("Triaged (Majors)", "RAT Room", 1),
    ("Triaged (Majors)", "Majors Bays", 2),
    ("Triaged (Majors)", "Majors Cubicles", 2),
    ("Triaged (Majors)", "Majors Corridor", 3),
    ("Triaged (Majors)", "Ambulance", 4),
    ("Triaged (Ambulatory)", "Ambulatory Triage", 1),
    ("Triaged (Ambulatory)", "Ambulatory Cubicles", 2),
    ("Triaged (Minors)", "Minors Triage", 1),
    ("Clerked (Resus)", "Resus Bays", 1),
    ("Clerked (Resus)", "Ambulance", 3),
    ("Clerked (Majors)", "Majors Bays", 1),
    ("Clerked (Majors)", "Majors Cubicles", 1),
    ("Clerked (Majors)", "Majors Corridor", 2),
    ("Clerked (Majors)", "Ambulance", 3),
    ("Clerked (Ambulatory)", "Ambulatory Cubicles", 1),
    ("Clerked (Minors)", "Treatment Rooms", 1),
    ("Decision to Admit (Resus)", "Resus Bays", 1),
    ("Decision to Admit (Resus)", "Ambulance", 3),
    ("Decision to Admit (Majors)", "Majors Bays", 1),
    ("Decision to Admit (Majors)", "Majors Cubicles", 1),
    ("Decision to Admit (Majors)", "Majors Corridor", 2),
    ("Decision to Admit (Majors)", "Ambulance", 3),
    ("Decision to Admit (Ambulatory)", "Ambulatory Cubicles", 1),
    ("Decision to Admit (Minors)", "Treatment Rooms", 1),
    ("Specialty Reviewed (Resus)", "Resus Bays", 1),
    ("Specialty Reviewed (Resus)", "Ambulance", 3),
    ("Specialty Reviewed (Majors)", "Majors Bays", 1),
    ("Specialty Reviewed (Majors)", "Majors Cubicles", 1),
    ("Specialty Reviewed (Majors)", "Majors Corridor", 2),
    ("Specialty Reviewed (Majors)", "Ambulance", 3),
    ("Specialty Reviewed (Ambulatory)", "Ambulatory Cubicles", 1),
    ("Specialty Reviewed (Minors)", "Treatment Rooms", 1),
    ("Senior Reviewed (Resus)", "Resus Bays", 1),
    ("Senior Reviewed (Resus)", "Ambulance", 3),
    ("Senior Reviewed (Majors)", "Majors Bays", 1),
    ("Senior Reviewed (Majors)", "Majors Cubicles", 1),
    ("Senior Reviewed (Majors)", "Majors Corridor", 2),
    ("Senior Reviewed (Majors)", "Ambulance", 3),
    ("Senior Reviewed (Ambulatory)", "Ambulatory Cubicles", 1),
    ("Senior Reviewed (Minors)", "Treatment Rooms", 1),
    ("Seen By Clinician/Treated (Resus)", "Resus Bays", 1),
    ("Seen By Clinician/Treated (Resus)", "Ambulance", 3),
    ("Seen By Clinician/Treated (Majors)", "Majors Bays", 1),
    ("Seen By Clinician/Treated (Majors)", "Majors Cubicles", 1),
    ("Seen By Clinician/Treated (Majors)", "Majors Corridor", 2),
    ("Seen By Clinician/Treated (Majors)", "Ambulance", 3),
    ("Seen By Clinician/Treated (Ambulatory)", "Ambulatory Cubicles", 1),
    ("Seen By Clinician/Treated (Minors)", "Treatment Rooms", 1)]

#######################LOCATION CAPACITY DATA#######################
location_capacity_data = {location[1]: "Unlimited"
                          for location in process_location_data}
location_capacity_data.update({"Resus Bays": 7,
                               "RAT Room": 2,
                               "Majors Bays": 10,
                               "Majors Cubicles": 9,
                               "Majors Corridor": 6,
                               "Treatment Rooms": 5,
                               "Minors Triage": 2,
                               "Ambulatory Triage": 3,
                               "Ambulatory Cubicles": 7,
                               "Reception": 2})
