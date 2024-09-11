from pathlib import Path
import data_cleaning_and_transformation as cleaning
import pathway_definitions as pathways
import process_durations as durations
import config
import pandas as pd
import numpy as np

if __name__ == "__main__":

    input_path = Path(r"./Events Data")
    output_path = Path(r"./Outputs/Baseline")

################################################################################
#-------------------------Read in and clense raw data--------------------------#
################################################################################

    # ---------------------- Events data
    events_raw = pd.read_csv(f"{input_path}/{'FN_NursingAssessmentsData.csv'}")
    events_quality = cleaning.drop_duplicates_and_anomaly_times_events_data(
                     events_raw, config.repeat_time_threshold,
                     config.remove_duplicate_staffid,
                     config.remove_duplicate_location)

    # ---------------------- Diagnostics data
    diagnostics_quality = None
    if config.include_diag_data:
        diagnostics_raw = pd.read_csv(f"{input_path}/{"FN_Diagnostics.csv"}")
        diagnostics_quality = (cleaning
                               .rename_columns_and_collapse_data_diagnostics(
                                   config.collapse_diagnostics_within_time,
                                   diagnostics_raw))

    # ---------------------- Obs data
    obs_quality = None
    if config.include_obs_data:
        obs_raw = pd.read_csv(f"{input_path}/{"FN_Obs.csv"}")
        if obs_raw is not None:
            obs_quality = cleaning.rename_columns_and_change_events_to_obs(
                          obs_raw)

    # --------------------- Admission data
    adm_status_raw = None
    if config.include_admission_data:
        adm_status_raw = pd.read_csv(f"{input_path}/{"FN_AdmissionStatus.csv"}")

    # ---------------------- Combine and Clense data
    events_quality = cleaning.main_cleanse_and_transform_data(
                     events_quality, adm_status_raw, obs_quality,
                     diagnostics_quality,
                     config.event_names_to_exclude_for_repetition,
                     config.excluded_event_names, config.locations_to_drop,
                     config.natural_order_for_processes,
                     config.include_spawn_end_events,
                     config.locations_pathway_map, config.keep_last_location,
                     config.admitted_map)

################################################################################
#-----------------------------Pathway Definitions------------------------------#
################################################################################

    # ---------------------- Calculate the number of patients making each
    #                        transition.
    transitions = pathways.add_reset_transitions(events_quality)

    # ----------------------- Definition Pathways generation for all data
    analysis_name = "Split by Pathway - All"
    filepath = output_path / "Pathways" / analysis_name
    pathways.main_generate_dfg_and_pathway_definitions(analysis_name,
        transitions, filepath, config.include_spawn_end_events,
        config.obs_splits, config.export_event_log_csv,
        config.export_log_to_csv_after_using_log_converter,
        config.pathways_wait_in_place, np.nan, False, False, "EventName",
        "Pathway")
    
    # ------------------------ Definition pathways generation with exclusion
    # ------------------------ under a certain threshold
    for threshold in range(1, 4):

        # ------------------------ Removed patient from data if they have a
        # ------------------------ transition below the threshold percentage
        analysis_name = f"Visits Exclusion {threshold}%"
        filepath = output_path / "Pathways" / analysis_name
        pathways.main_generate_dfg_and_pathway_definitions(analysis_name,
                 transitions, filepath, config.include_spawn_end_events,
                 config.obs_splits, config.export_event_log_csv,
                 config.export_log_to_csv_after_using_log_converter,
                 config.pathways_wait_in_place, threshold, True, False,
                 "EventName", "Pathway")
        
        # ------------------------ Remove transitions that fall below the
        # ------------------------ threshold percentage
        directory_path = f"Removed transitions below {threshold}%"
        filepath = output_path / "Pathways" / directory_path
        pathways.main_generate_dfg_and_pathway_definitions(directory_path,
                 transitions, filepath, config.include_spawn_end_events,
                 config.obs_splits, config.export_event_log_csv,
                 config.export_log_to_csv_after_using_log_converter,
                 config.pathways_wait_in_place, threshold, False, True,
                 "EventName", "Pathway")

################################################################################
#------------------------------Process Durations-------------------------------#
################################################################################

    # ------------------------ Calculate the duration of each event in the data
    event_diffs = durations.add_difference_in_minutes_to_durations(
                  events_quality, config.where_duration_should_be_0)
    
    # ------------------------ Generate and output histograms and process
    # ------------------------ durations for different filterings/scenarios
    analysis_name = "Max threshold 2 hours and including 100 percentile"
    durations.main_generate_histogram_and_process_durations(analysis_name,
              event_diffs, "Event (Pathway)", output_path, config.plots,
              [durations.within_threshold_diff(120)])
    
    # ------------------------
    analysis_name = "Max threshold 2 hours, 100 perc, between 8am and 10pm"
    durations.main_generate_histogram_and_process_durations(analysis_name,
              event_diffs, "Event (Pathway)", output_path, config.plots,
              [durations.within_threshold_diff(120),
               durations.only_daytime_events])
    
    # ------------------------
    analysis_name = "Max threshold 14 hours,"\
                    f" and {int(config.quantile_threshold*100)} percentile"
    durations.main_generate_histogram_and_process_durations(analysis_name,
              event_diffs, "Event (Pathway)", output_path, config.plots,
              [durations.within_threshold_diff(840),
               durations.within_diff_quantile(config.quantile_threshold)])
    
    # ------------------------
    analysis_name = "Max threshold 14 hours, "\
                    f"{int(config.quantile_threshold*100)} perc,"\
                    " between 8am and 10pm"
    durations.main_generate_histogram_and_process_durations(analysis_name,
              event_diffs, "Event (Pathway)", output_path, config.plots,
              [durations.within_threshold_diff(840),
               durations.within_diff_quantile(config.quantile_threshold),
               durations.only_daytime_events])
