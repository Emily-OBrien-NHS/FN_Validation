from utils import sort_events
from config import WALK_IN, SPAWN, REMOVED, WAITING_FOR_BED
import pandas as pd

#-----------------------------------initial cleaning functions on raw data
def drop_duplicates_and_anomaly_times_events_data(events_raw,
                                                  repeat_time_threshold,
                                                  remove_duplicate_staffid,
                                                  remove_duplicate_location):
    """
    Args:
        events_raw (pd.DataFrame): raw events data frame.
        repeat_time_threshold (int): the number of minutes that should be used
        to filter out repeated events (e.g. two triage on the same patient
        within n minutes is likely a duplicate and should be discounted).
        remove_duplicate_staffid (bool): Flag to consider if duplicate events
        (within n minutes) should be removed if done by the same staff member.
        remove_duplicate_location (bool): Flag to consider if duplicate events
        (within n minutes) should be removed if done in the same location.

    Returns:
        pd.DataFrame: version of events df with duplicates removed.
    """

    #Remove fully duplicated entries, format datetime column and remove
    #any data from before 2018-04-01
    events_quality = events_raw.copy().drop_duplicates()
    events_quality["EventTime"] = pd.to_datetime(events_quality["EventTime"],
                                                 format="%d/%m/%Y %H:%M")
    events_quality = events_quality.dropna(subset=["EventTime"])

    events_quality = events_quality.loc[~(events_quality["EventTime"]
                                        < pd.to_datetime("2018-04-01"))].copy()

    #Add diff column of the time between duplicate events for the same VisitId.
    #Fill in NaT with the threshold+1 so these don't get filtered out.
    events_quality = events_quality.sort_values(by=['VisitId', 'EventName',
                                                    'EventTime'])
    events_quality['diff'] = ((events_quality.groupby(['VisitId', 'EventName'])
                              ['EventTime'].diff().infer_objects(copy=False)
                              .fillna(pd.Timedelta('NaT')).dt.seconds / 60)
                              .fillna(repeat_time_threshold+1))
    #Empty list of flag cols, populate if any of the flags are True
    flagcols = []
    if remove_duplicate_staffid:
        events_quality['DuplicateStaffId'] = (events_quality['EventStaffId']
                                    == events_quality['EventStaffId'].shift(1))
        flagcols.append('DuplicateStaffId')

    if remove_duplicate_location:
        events_quality['DuplicateLocation'] = (events_quality['EventLocation']
                                    == events_quality['EventLocation'].shift(1))
        flagcols.append('DuplicateLocation')

    #if there are flags set to true, add the filter for these into the mask to
    #drop duplicate events
    if flagcols:
        mask = ~((events_quality['diff'] < repeat_time_threshold)
                 & (events_quality[flagcols].all(axis=1)))
    else:
        mask = ~(events_quality['diff'] < repeat_time_threshold)

    #use mask to remove duplicate events within n minutes of each other
    events_quality = events_quality.loc[mask, ['VisitId', 'EventName',
                                               'EventTime', 'EventStaffId',
                                               'EventLocation']].copy()
    return events_quality


def rename_columns_and_collapse_data_diagnostics(
        collapse_diagnostics_rows_within_time_of, diagnostics_raw):
    """
    Args:
        collapse_diagnostics_rows_within_time_of (pd.Timedelta): time between
        diagnostic event times to remove duplicates within n minutes.
        diagnostics_raw (pd.DataFrame): raw diagnostics data frame.

    Returns:
        pd.DataFrame: version of diagnostics with duplicates removed.
    """

    diagnostics_quality = diagnostics_raw[["VisitID", "Request DateTime",
                                           "ItemMasterCategory"]].copy()
    diagnostics_quality = diagnostics_quality.rename(columns={
                                            "ItemMasterCategory": "EventName",
                                            "Request DateTime": "EventTime",
                                            "VisitID": "VisitId"})
    diagnostics_quality["EventTime"] = pd.to_datetime(
                    diagnostics_quality["EventTime"], format="%d/%m/%Y %H:%M")

    # collapse similar rows within a short time for a patient to one row
    mask = (diagnostics_quality.sort_values(by=["VisitId", "EventTime"])
            .groupby("VisitId")["EventTime"].diff()
            .lt(collapse_diagnostics_rows_within_time_of).sort_index())

    # index the df with this mask to remove the duplicates
    diagnostics_quality = diagnostics_quality.loc[~mask].copy()
    return diagnostics_quality


def rename_columns_and_change_events_to_obs(obs_raw):
    """
    Args:
        obs_raw (pd.DataFrame): raw obs data frame.

    Returns:
       pd.DataFrame: a cleaned version of the obs data frame.
    """
    obs_quality = obs_raw.copy()
    obs_quality = obs_raw.rename(columns={"ChartType": "EventName",
                                          "ChartDateTime": "EventTime",
                                          "VisitID": "VisitId"})
    obs_quality["EventTime"] = pd.to_datetime(obs_quality["EventTime"],
                                              format="%d/%m/%Y %H:%M")
    obs_quality["EventName"] = "Observations"
    return obs_quality

#------------------------Cleaning functions used in main cleaning function
def remove_repeats_of_events_that_should_not_be_repeated(events_quality,
                                        event_name_to_exclude_for_repetition):
    """
    Args:
        events_quality (pd.DataFrame): clensed events data frame.
        event_name_to_exclude_for_repetition (list[str]): list of events that
        should not repeat.

    Returns:
        pd.DataFrame: clensed events data frame with duplicates of non repeat
        events removed.
    """
    mask_for_dropping_duplicate_events = (events_quality["EventName"]
                                    .isin(event_name_to_exclude_for_repetition)
                                    & (events_quality.duplicated(subset=
                                                                 ["VisitId",
                                                                  "EventName"],
                                                                keep="first")))
    events_quality = (events_quality
                      .loc[~mask_for_dropping_duplicate_events].copy())
    return events_quality


def remove_events_after_discharged(events_quality):
    """
    Args:
        events_quality (pd.DataFrame): clensed events data frame.

    Returns:
        pd.DataFrame: events data frame with any events after discharge removed.
    """
    events_quality = events_quality.sort_values(by=["VisitId", "EventTime"])
    target_event = "Discharged"
    events_quality["target_event"] = events_quality["EventName"] == target_event
    events_quality["after_target_event"] = (events_quality.groupby("VisitId")
                                        ["target_event"].cumsum().astype(bool))
    events_quality = (events_quality.loc[~(events_quality["after_target_event"])
                                         | (events_quality["target_event"])]
                                         .copy().drop(columns=["target_event",
                                                        "after_target_event"]))
    return events_quality


def augmenting_admittance_data(adm_status_raw, events_quality):

    """_summary_
    Args:
        adm_status_raw (pd.DataFrame): admitted status raw data frame.
        events_quality (pd.DataFrame): clensed events data frame.
    Returns:
        pd.DataFrame: clensed events data frame.
    """
    adm_status_quality = adm_status_raw.rename(columns={"AttendanceID"
                                                        : "VisitId"})
    adm_status_quality = adm_status_quality.loc[adm_status_quality["Adm"]
                                                != "Non-Admitted",
                                                ["VisitId", "Adm"]].copy()
    events_quality = events_quality.merge(adm_status_quality, on="VisitId",
                                          how="left")
    events_quality.loc[(events_quality["EventName"] == "Discharged")
                       & (~pd.isnull(events_quality["Adm"])),
                       "EventName"] = events_quality["Adm"]
    events_quality = events_quality.drop(["Adm"], axis=1)
    return events_quality


def merge_data(events_quality, diagnostics_quality, obs_quality):
    """
    Args:
        events_quality (pd.DataFrame): clensed events dataframe
        diagnostics_quality (Optional[pd.DataFrame], optional): clensed
        diagnostics data to merge onto events data. Defaults to None.
        obs_quality (Optional[pd.DataFrame], optional): clensed obs data to
        merge onto events data. Defaults to None.

    Returns:
        pd.DataFrame: clensed events data merged onto other required data
    """
    if obs_quality is not None:
        events_quality = pd.concat([events_quality, obs_quality])
    if diagnostics_quality is not None:
        events_quality = pd.concat([events_quality, diagnostics_quality])
    events_quality["EventTime"] = pd.to_datetime(events_quality["EventTime"],
                                                  format="%d/%m/%Y %H:%M")
    return events_quality


def remove_excluded_events_and_locations(events_quality, excluded_event_names,
                                         locations_to_drop):
    """
    Args:
        events_quality (pd.DataFrame): clensed events dataframe.
        excluded_event_names (list[str]): list of event names to exclude.
        locations_to_drop (list[str]): list of locations to exclude.

    Returns:
        pd.DataFrame: clensed events dataframe.
    """
    #remove excluded event names
    events_quality = events_quality.loc[~events_quality["EventName"]
                                        .isin(excluded_event_names)].copy()
    #remove visit ids with excluded event locations
    mask = events_quality.loc[events_quality["EventLocation"]
                              .isin(locations_to_drop)
                              | events_quality["EventLocation"].isna()].copy()
    visit_ids_to_exclude = mask["VisitId"].unique()
    events_quality = events_quality.loc[~events_quality["VisitId"]
                                        .isin(visit_ids_to_exclude)].copy()

    return events_quality


def set_location_for_ambulance_arrival(events_quality):
    """
    Args:
        events_quality (pd.DataFrame): clensed events dataframe.

    Returns:
        pd.DataFrame: clensed events data frame with Ambulance as event location
        where required.
    """
    events_quality.loc[(pd.isnull(events_quality["EventLocation"]))
                       & (events_quality["EventName"] == "Ambulance Arrival"),
                       "EventLocation"] = "Ambulance"

    return events_quality


def forward_fill_on_locations(events_quality):
    """
    Args:
        events_quality (pd.DataFrame): clensed events dataframe.
    Returns:
        pd.DataFrame: clensed events dataframe with missing event locations
        filled in.
    """
    events_quality["EventLocation"] = (events_quality.groupby("VisitId")
                                       ["EventLocation"].ffill())

    return events_quality


def mapping_of_natural_order(events_quality, natural_order_for_processes):
    """
    Args:
        events_quality (pd.DataFrame): clensed events dataframe.
        natural_order_for_processes (dict[str, int]): dictionary of the ideal
        order of processes.

    Returns:
        pd.DataFrame: clensed dataframe with an order number assigned to events.
    """

    events_quality["natural_order"] = (events_quality["EventName"]
                                       .map(natural_order_for_processes))

    return events_quality


def add_walk_in_for_non_ambulance_arrivals(events_quality):
    """
    Args:
        events_quality (pd.DataFrame): clensed events dataframe.

    Returns:
        pd.DataFrame: clensed events dataframe with walkin for non ambulance
        arrivals.
    """
    # groupby to get first event for each patient
    walk_in_filter = (events_quality.groupby("VisitId", as_index=False)
                      [["EventTime", "EventName"]].first())

    # filter out all ambulance arrivals
    walk_in_filter = walk_in_filter.loc[walk_in_filter["EventName"]
                                        != "Ambulance Arrival"].copy()

    # groupby to get min values just for walk-ins
    only_walk_ins_first_event = (walk_in_filter
                                 .groupby("VisitId",as_index=False)
                                 [["EventTime"]].min())

    # copy df and add the walk in event
    only_walk_ins_first_event["EventName"] = WALK_IN
    # add back to main df
    events_quality = pd.concat([events_quality, only_walk_ins_first_event])

    return events_quality


def add_wait_for_beds_for_admitted_patients(events_quality):
    """
    Args:
        events_quality (pd.DataFrame): clensed events dataframe.

    Returns:
        pd.DataFrame: clensed events dataframe with wait for bed as an event for
        admitted patients.
    """

    events_quality = sort_events(events_quality)
    #Get last event for each patient
    lastevent_filter_groupby = (events_quality
                                 .groupby("VisitId", as_index=False)
                                 [["EventTime", "EventName"]].last())
    #remove any discharged patients
    admitted_patients = lastevent_filter_groupby.loc[
                                 lastevent_filter_groupby["EventName"]
                                 != "Discharged"].copy()
    #Add in the extra wait for bed event
    admitted_patients["EventName"] = (WAITING_FOR_BED + " - ("
                                      + admitted_patients["EventName"] + ")")

    events_quality = pd.concat([events_quality, admitted_patients])

    return events_quality


def add_spawn_end_events(events_quality, natural_order_for_processes):
    """
    Args:
        events_quality (pd.DataFrame): clensed events dataframe.

    Returns:
        pd.DataFrame: clensed events dataframe with spawned and removed events.
    """
    # groupby to get min values - spawn
    min_df = events_quality.groupby("VisitId", as_index=False)["EventTime"].first()
    min_df["EventName"] = SPAWN
    min_df["natural_order"] = natural_order_for_processes[SPAWN]

    # groupby to get max values - end
    max_df = events_quality.groupby("VisitId", as_index=False)["EventTime"].max()
    max_df["EventName"] = REMOVED
    max_df["natural_order"] = natural_order_for_processes[REMOVED]

    # concat these to main dataframe
    events_quality = pd.concat([events_quality, min_df, max_df])
    events_quality = sort_events(events_quality)
    events_quality["EventLocation"] = events_quality["EventLocation"].ffill()

    return events_quality


def keep_last_location_of_patient(events_quality):
    """
    Args:
        events_quality (pd.DataFrame): clensed events dataframe.
    Returns:
        pd.DataFrame: clensed events dataframe with last event location kept.
    """
    events_quality["EventLocation"] = (events_quality.groupby("VisitId")
                                       ["EventLocation"].transform("last"))

    return events_quality


def map_locations_to_triage_category_and_create_pathway_column(events_quality,
                                                         locations_pathway_map):
    """
    Args:
        events_quality (pd.DataFrame): clensed events quality dataframe.
        locations_pathway_map (dict[str, str]): dictionary to map locations to
        their pathway.
        natural_order_for_processes (dict[str, int]): dictionary to map 
        processes to their order.

    Returns:
        pd.DataFrame: clensed events quality dataframe woth pathway addded.
    """
    #Map the pathways to the event locations
    events_quality["Pathway"] = (events_quality["EventLocation"]
                                 .map(locations_pathway_map))
    #Create a column of the event name and pathway
    events_quality["Event (Pathway)"] = (events_quality["EventName"]
                                       + " (" + events_quality["Pathway"] + ")")
    events_quality = sort_events(events_quality)

    return events_quality
