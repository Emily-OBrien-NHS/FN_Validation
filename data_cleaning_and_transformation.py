import pandas as pd

################################################################################
#--------------------initial cleaning functions on raw data--------------------#
################################################################################

def drop_duplicates_and_anomaly_times_events_data(events_raw,
                                                  repeat_time_threshold,
                                                  remove_duplicate_staffid,
                                                  remove_duplicate_location):
    #Function to remove duplicates and anomoly time events data.
    #----
    #Remove fully duplicated entries, format datetime column and remove
    #any data from before 2018-04-01 and missing data.
    events_quality = events_raw.copy().drop_duplicates()
    events_quality["EventTime"] = pd.to_datetime(events_quality["EventTime"],
                                                 format="%d/%m/%Y %H:%M")
    events_quality = events_quality.dropna(subset=["EventTime"])
    events_quality = events_quality.loc[(events_quality["EventTime"]
                                        >= pd.to_datetime("2018-04-01"))].copy()
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
        flag_col = 'DuplicateStaffId'
        events_quality[flag_col] = (events_quality['EventStaffId']
                                    == events_quality['EventStaffId'].shift(1))
        flagcols.append(flag_col)
    if remove_duplicate_location:
        flag_col = 'DuplicateLocation'
        events_quality[flag_col] = (events_quality['EventLocation']
                                    == events_quality['EventLocation'].shift(1))
        flagcols.append(flag_col)
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
    #Function to rename columns and collapse similar events in a short time
    #frame to one event for diagnostics data.
    #----
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
    #Function to format obs data
    #----
    obs_quality = obs_raw.copy()
    obs_quality = obs_raw.rename(columns={"ChartType": "EventName",
                                          "ChartDateTime": "EventTime",
                                          "VisitID": "VisitId"})
    obs_quality["EventTime"] = pd.to_datetime(obs_quality["EventTime"],
                                              format="%d/%m/%Y %H:%M")
    obs_quality["EventName"] = "Observations"
    return obs_quality

################################################################################
#---------------Cleaning functions used in main cleaning function--------------#
################################################################################

def sort_events(events):
    #Function to sort clensed events dataframe (used several times)
    events = events.sort_values(by=["VisitId", "EventTime", "natural_order"])
    return events

def remove_repeats_of_events_that_should_not_be_repeated(events_quality,
                                        event_names_to_exclude_for_repetition):
    #Function to remove any repeats of events that should not be repeated
    duplicate_events_mask = (events_quality["EventName"]
                             .isin(event_names_to_exclude_for_repetition)
                            & (events_quality
                               .duplicated(subset=["VisitId", "EventName"],
                                           keep="first")))
    events_quality = events_quality.loc[~duplicate_events_mask].copy()
    return events_quality

def remove_events_after_discharged(events_quality):
    #Remove any events that happen after discharge
    events_quality = events_quality.sort_values(by=["VisitId", "EventTime"])
    events_quality["target_event"] = events_quality["EventName"] == "Discharged"
    events_quality["after_target_event"] = (events_quality.groupby("VisitId")
                                        ["target_event"].cumsum().astype(bool))
    events_quality = (events_quality
                      .loc[~(events_quality["after_target_event"])
                           | (events_quality["target_event"])].copy()
                      .drop(columns=["target_event", "after_target_event"]))
    return events_quality

def remove_senior_review_if_seen_by(events_quality):
    #function to remove senior review events if seen by event at the same time/
    #staff id/ location.
    #----
    #The below works because alphabetically, senior reviewed is after seen by.
    rem = (events_quality
           .loc[events_quality['EventName'].isin(['Seen By Clinician/Treated',
                                                  'Senior Reviewed'])]
           .sort_values(by=['EventName'])
           .duplicated(subset=['VisitId', 'EventTime', 'EventStaffId',
                               'EventLocation'], keep='first'))
    #Remove these events
    events_quality = (events_quality.loc[~events_quality.index
                                         .isin(rem[rem].index)].copy())
    return events_quality

def augmenting_admittance_data(adm_status_raw, events_quality):
    #Function to add in admitted data if provided.
    adm_status_quality = adm_status_raw.rename(columns={"AttendanceID"
                                                        : "VisitId"})
    adm_status_quality = (adm_status_quality
                          .loc[adm_status_quality["Adm"] != "Non-Admitted",
                          ["VisitId", "Adm"]].copy())
    events_quality = events_quality.merge(adm_status_quality, on="VisitId",
                                          how="left")
    events_quality.loc[(events_quality["EventName"] == "Discharged")
                       & (~pd.isnull(events_quality["Adm"])),
                       "EventName"] = events_quality["Adm"]
    events_quality = events_quality.drop(["Adm"], axis=1)
    return events_quality

def merge_data(events_quality, diagnostics_quality, obs_quality):
    #Function to add in obs and diagnostics data if the data is provided.
    if obs_quality is not None:
        events_quality = pd.concat([events_quality, obs_quality])
    if diagnostics_quality is not None:
        events_quality = pd.concat([events_quality, diagnostics_quality])
    events_quality["EventTime"] = pd.to_datetime(events_quality["EventTime"],
                                                  format="%d/%m/%Y %H:%M")
    return events_quality

def set_location_for_ambulance_arrival(events_quality):
    #Function to add Ambulance location for Ambulance Arrivals
    events_quality.loc[(pd.isnull(events_quality["EventLocation"]))
                       & (events_quality["EventName"] == "Ambulance Arrival"),
                       "EventLocation"] = "Ambulance"
    return events_quality

def forward_fill_on_locations(events_quality):
    #Function to forward fill missing event locations
    events_quality["EventLocation"] = (events_quality.groupby("VisitId")
                                       ["EventLocation"].ffill())
    return events_quality

def remove_excluded_events_and_locations(events_quality, excluded_event_names,
                                         locations_to_drop):
    #Function to remove excluded events and locations
    #----
    #remove excluded event names
    events_quality = events_quality.loc[~events_quality["EventName"]
                                        .isin(excluded_event_names)].copy()
    #remove visit ids with excluded event locations
    visit_ids_to_exclude = (events_quality
                            .loc[events_quality["EventLocation"]
                                 .isin(locations_to_drop)
                                | events_quality["EventLocation"].isna(),
                                'VisitId'].unique())
    events_quality = events_quality.loc[~events_quality["VisitId"]
                                        .isin(visit_ids_to_exclude)].copy()
    return events_quality

def mapping_of_natural_order(events_quality, natural_order_for_processes):
    #Function to map events to their order ready for sorting.
    events_quality["natural_order"] = (events_quality["EventName"]
                                       .map(natural_order_for_processes))
    return events_quality

def add_walk_in_for_non_ambulance_arrivals(events_quality):
    #Function to add in walkin arrivals for any non-ambulance arrival.
    #----
    # groupby to get first event for each patient
    walk_in_filter = (events_quality.groupby("VisitId", as_index=False)
                      [["EventTime", "EventName"]].first())
    # filter out all ambulance arrivals
    walk_in_filter = walk_in_filter.loc[walk_in_filter["EventName"]
                                        != "Ambulance Arrival"].copy()
    # groupby to get min values just for walk-ins
    only_walk_ins_first_event = (walk_in_filter
                                 .groupby("VisitId", as_index=False)
                                 ["EventTime"].min())
    # copy df and add the walk in event
    only_walk_ins_first_event["EventName"] = "Walk-In"
    # add back to main df
    events_quality = pd.concat([events_quality, only_walk_ins_first_event])
    return events_quality


def add_wait_for_beds_for_admitted_patients(events_quality, admitted_map):
    #Add wait for bed events if there is an event after discharge
    #----
    #Get last event for each patient
    last_event = (events_quality.groupby("VisitId", as_index=False)
                  [["EventTime", "EventName"]].last())
    #remove any discharged patients
    admitted_patients = (last_event.loc[last_event["EventName"] != "Discharged"]
                         .copy())
    #Add in the extra wait for bed event
    admitted_patients["EventName"] = ("Wait for Bed - "
                                      + admitted_patients["EventName"])
    #rename admitted events to just admitted
    events_quality['EventName'] = (events_quality['EventName']
                                   .replace(admitted_map))
    #Concat the admitted events to the events dataframe
    events_quality = pd.concat([events_quality, admitted_patients])
    return events_quality

def add_spawn_end_events(events_quality, natural_order_for_processes):
    #Function to add spawn and end events
    #----
    # groupby to get min values - spawn
    min_df = (events_quality.groupby("VisitId", as_index=False)
              ["EventTime"].first())
    min_df["EventName"] = "Spawn"
    min_df["natural_order"] = natural_order_for_processes["Spawn"]
    # groupby to get max values - end
    max_df = (events_quality.groupby("VisitId", as_index=False)
              ["EventTime"].max())
    max_df["EventName"] = "Removed"
    max_df["natural_order"] = natural_order_for_processes["Removed"]
    # concat these to main dataframe
    events_quality = pd.concat([events_quality, min_df, max_df])
    events_quality = sort_events(events_quality)
    events_quality["EventLocation"] = events_quality["EventLocation"].ffill()
    return events_quality

def keep_last_location_of_patient(events_quality):
    #Function to keep the last location of each patient as the location
    events_quality["EventLocation"] = (events_quality.groupby("VisitId")
                                       ["EventLocation"].transform("last"))
    return events_quality

def map_locations_and_create_pathway_column(events_quality,
                                                         locations_pathway_map):
    #Function to map event locations to a pathway and create the pathway column
    #of EventName (Pathway)
    #----
    #Map the pathways to the event locations
    events_quality["Pathway"] = (events_quality["EventLocation"]
                                 .map(locations_pathway_map))
    #Create a column of the event name and pathway
    events_quality["Event (Pathway)"] = (events_quality["EventName"]
                                       + " (" + events_quality["Pathway"] + ")")
    return events_quality

################################################################################
#----------------------------main cleaning function----------------------------#
################################################################################

def main_cleanse_and_transform_data(events_quality, adm_status_raw, obs_quality,
                                    diagnostics_quality,
                                    event_names_to_exclude_for_repetition,
                                    excluded_event_names, locations_to_drop,
                                    natural_order_for_processes,
                                    include_spawn_end_events,
                                    locations_pathway_map, keep_last_location,
                                    admitted_map):
    #Main function to apply all cleaning and transformation functions onto the
    #input data
    #----
    #remove repeated or 'non real' events
    events_quality = (remove_repeats_of_events_that_should_not_be_repeated(
                       events_quality, event_names_to_exclude_for_repetition))
    events_quality = remove_events_after_discharged(events_quality)
    events_quality = remove_senior_review_if_seen_by(events_quality)
    #Add admissions, diagnosis and obs data if required
    if adm_status_raw is not None:
        events_quality = augmenting_admittance_data(adm_status_raw,
                                                    events_quality)
    events_quality = merge_data(events_quality, diagnostics_quality,
                                         obs_quality)
    #fill in missing locations and remove locations/enets for removal
    events_quality = set_location_for_ambulance_arrival(events_quality)
    events_quality = forward_fill_on_locations(events_quality)
    events_quality = remove_excluded_events_and_locations(events_quality,
                                                          excluded_event_names,
                                                          locations_to_drop)
    #Sort events data and add in walkin and wait for bed events   
    events_quality = mapping_of_natural_order(events_quality,
                                              natural_order_for_processes)
    events_quality = sort_events(events_quality)
    events_quality = add_walk_in_for_non_ambulance_arrivals(events_quality)
    events_quality = sort_events(events_quality)
    events_quality = add_wait_for_beds_for_admitted_patients(events_quality,
                                                             admitted_map)
    events_quality = mapping_of_natural_order(events_quality,
                                              natural_order_for_processes)
    #If true, add in spawn and end events
    if include_spawn_end_events:
        events_quality = add_spawn_end_events(events_quality,
                                              natural_order_for_processes)
    #If true, replace all locations with patient's final location
    if keep_last_location:
        events_quality = keep_last_location_of_patient(events_quality)
    #Add pathway column and sort events
    events_quality = map_locations_and_create_pathway_column(events_quality,
                                                        locations_pathway_map)
    events_quality = sort_events(events_quality)
    return events_quality
