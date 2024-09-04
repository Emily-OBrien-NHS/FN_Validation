from pathlib import Path
from utils import load_data
from event_filtering_functions import (
    exclude_patients_with_uncommon_transitions_below_threshold,
    within_diff_quantile,
    within_threshold_diff,
    only_daytime_events)
from main_data_cleaning_function import cleanse_and_transform_data
import data_cleaning_and_transformation as cleaning
import pathway_definitions as pathways
import process_durations as durations
import config
import pandas as pd

if __name__ == "__main__":

    output_path = Path(r"./Outputs")
    path_to_read_data = Path(r"./Events Data")
    # ---------------------- Read in and clense raw data
    # ---------------------- Events data
    #events_raw = load_data("FN_Events.csv")
    events_raw = load_data('FN_NursingAssessmentsData.csv')
    events_quality = cleaning.drop_duplicates_and_anomaly_times_events_data(
        events_raw, config.repeat_time_threshold,
        config.remove_duplicate_staffid, config.remove_duplicate_location)

    # ---------------------- Diagnostics data
    if config.include_diag_data:
        diagnostics_raw = load_data("FN_Diagnostics.csv")
        diagnostics_quality = cleaning.rename_columns_and_collapse_data_diagnostics(
            pd.Timedelta("5m"), diagnostics_raw)
    else:
        diagnostics_quality = None

    # ---------------------- Obs data
    obs_quality = None
    if config.include_obs_data:
        obs_raw = load_data("FN_Obs.csv")
        if obs_raw is not None:
            obs_quality = cleaning.rename_columns_and_change_events_to_obs(obs_raw)

    # --------------------- Admission data
    adm_status_raw = None
    if config.include_admission_data:
        adm_status_raw = load_data("FN_AdmissionStatus.csv")

    # ---------------------- Clense data
    events_quality = cleanse_and_transform_data(events_quality, adm_status_raw,
                                                obs_quality, diagnostics_quality,
                                                config.excluded_event_names,
                                                config.locations_to_drop,
                                                config.natural_order_for_processes,
                                                config.include_spawn_end_events,
                                                config.locations_pathway_map,
                                                config.keep_last_location)

    # ---------------------- Calculate the number of patients making each
    #                        transition.
    transitions = pathways.add_reset_transitions(events_quality)

    # ----------------------- Definition Pathways generation for all data
    # Define the directory_path name for outputs
    analysis_name = "Split by Pathway - All"
    filepath = output_path / "Pathways" / analysis_name

    pathways.generate_and_output_dfg_and_pathway_definition(analysis_name,
             transitions, filepath, config.include_spawn_end_events,
             config.obs_splits, config.export_event_log_csv,
             config.export_log_to_csv_after_using_log_converter,
             config.location_capacity_data, config.process_location_data,
             config.event_names_based_process_requirements,
             config.pathways_wait_in_place, "EventName", "Pathway")

    # ------------------------ Definition pathways generation with exclusion
    #                          under a certain threshold

    for threshold in range(2, 4):
        analysis_name = f"Visits Exclusion {threshold}%"
        filepath = output_path / "Pathways" / analysis_name

        pathways.generate_and_output_dfg_and_pathway_definition(analysis_name,
                 transitions, filepath, config.include_spawn_end_events,
                 config.obs_splits, config.export_event_log_csv,
                 config.export_log_to_csv_after_using_log_converter,
                 config.location_capacity_data, config.process_location_data,
                 config.event_names_based_process_requirements,
                 config.pathways_wait_in_place, "EventName", "Pathway",
                 [exclude_patients_with_uncommon_transitions_below_threshold(threshold)])

    for threshold in range(1, 4):
        directory_path = f"Removed transitions below {threshold}%"
        filepath = output_path / "Pathways" / directory_path

        pathways.generate_and_output_dfg_and_pathway_definition(directory_path,
                 transitions, filepath, config.include_spawn_end_events,
                 config.obs_splits, config.export_event_log_csv,
                 config.export_log_to_csv_after_using_log_converter,
                 config.location_capacity_data, config.process_location_data,
                 config.event_names_based_process_requirements,
                 config.pathways_wait_in_place, "EventName", "Pathway",
                 [],
                 pathways.remove_transitions_below_percentage_in_pathway_definitions(threshold))   


    # ------------------------------------- Process Durations
    event_diffs = durations.add_difference_in_minutes_to_durations(events_quality)

    analysis_name = "Max threshold 2 hours and including 100 percentile"
    durations.generate_and_output_histogram_and_process_durations(analysis_name,
              event_diffs, "Event (Pathway)", output_path, config.plots,
              [within_threshold_diff(120)])

    analysis_name = "Max threshold 14 hours, and 97 percentile"
    durations.generate_and_output_histogram_and_process_durations(analysis_name,
              event_diffs, "Event (Pathway)", output_path, config.plots,
              [within_threshold_diff(840), within_diff_quantile(0.97)])

    analysis_name = "Max threshold 14 hours, 97 perc, between 8am and 10pm"
    durations.generate_and_output_histogram_and_process_durations(analysis_name,
              event_diffs, "Event (Pathway)", output_path, config.plots,
              [within_threshold_diff(840), within_diff_quantile(0.97),
               only_daytime_events])

    analysis_name = "Max threshold 2 hours, 100 perc, between 8am and 10pm"
    durations.generate_and_output_histogram_and_process_durations(analysis_name,
              event_diffs, "Event (Pathway)", output_path, config.plots,
              [within_threshold_diff(120), only_daytime_events])

