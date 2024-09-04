from utils import sort_events
from config import event_names_to_exclude_for_repetition
import data_cleaning_and_transformation as cleaning
#import pandas as pd


def cleanse_and_transform_data(events_quality, adm_status_raw, obs_quality,
                               diagnostics_quality, excluded_event_names,
                               locations_to_drop, natural_order_for_processes,
                               include_spawn_end_events, locations_pathway_map,
                               keep_last_location):
    """
    Args:
        events_quality (pd.DataFrame): clensed events data frame.
        adm_status_raw (Optional[pd.DataFrame]): raw admission status data frame.
        obs_quality (Optional[pd.DataFrame]): clensed obs dataframe.
        diagnostics_quality (Optional[pd.DataFrame]): clensed diagnostics
        dataframe.
        excluded_event_names (list[str]): list of event names to exclude.
        locations_to_drop (list[str]): list of locations to exclude
        natural_order_for_processes (dict[str, int]): dictionary of order of 
        processes.
        include_spawn_end_events (bool): flag to include or exclude spawn end
        events.
        locations_pathway_map (dict[str, str]): dictionary to map locations to
        their pathway.
        keep_last_location (bool): flag to include last location.

    Returns:
        pd.DataFrame: a clensed events dataframe after applying all of the
        clensing steps.
    """

    events_quality = cleaning.remove_repeats_of_events_that_should_not_be_repeated(
                     events_quality, event_names_to_exclude_for_repetition)
    
    events_quality = cleaning.remove_events_after_discharged(events_quality)

    events_quality = cleaning.remove_senior_review_if_seen_by(events_quality)
    
    if adm_status_raw is not None:
        events_quality = cleaning.augmenting_admittance_data(adm_status_raw,
                                                             events_quality)

    events_quality = cleaning.merge_data(events_quality, diagnostics_quality,
                                         obs_quality)

    events_quality = cleaning.set_location_for_ambulance_arrival(events_quality)

    events_quality = cleaning.forward_fill_on_locations(events_quality)

    events_quality = cleaning.remove_excluded_events_and_locations(
                     events_quality, excluded_event_names, locations_to_drop)

    events_quality = cleaning.mapping_of_natural_order(events_quality,
                                                    natural_order_for_processes)

    events_quality = sort_events(events_quality)

    events_quality = cleaning.add_walk_in_for_non_ambulance_arrivals(events_quality)

    events_quality = cleaning.add_wait_for_beds_for_admitted_patients(events_quality)

    events_quality = cleaning.mapping_of_natural_order(events_quality,
                                                       natural_order_for_processes)

    if include_spawn_end_events:
        events_quality = cleaning.add_spawn_end_events(events_quality,
                                                       natural_order_for_processes)

    if keep_last_location:
        events_quality = cleaning.keep_last_location_of_patient(events_quality)

    events_quality = (cleaning
                    .map_locations_to_triage_category_and_create_pathway_column(
                    events_quality, locations_pathway_map))

    return events_quality
