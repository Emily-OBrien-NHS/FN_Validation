from config import REMOVED, SPAWN
import pandas as pd


def only_daytime_events(dataframe):
    """
    Args:
        dataframe (pd.DataFrame): Dataframe of events.

    Returns:
        pd.DataFrame: dataframe of only daytime events.
    """
    eight_am_to_10_pm = dataframe.copy()
    eight_am_to_10_pm.index = dataframe["EventTime"]
    eight_am_to_10_pm = eight_am_to_10_pm.between_time("08:00:00", "20:00:00")
    result = dataframe.loc[dataframe["EventTime"]
                           .isin(eight_am_to_10_pm.index)].copy()

    return result


def within_threshold_diff(max_diff_minutes_for_durations):
    """
    Args:
        max_diff_minutes_for_durations (float): number of minutes to use to
        filter out the periods between events.
    """
    def remove_data_under_certain_hours(dataframe):
        """
        Args:
            dataframe (pd.DataFrame): events dataframe

        Returns:
            pd.DataFrame: events dataframe with differences between events
            over the threshold removed.
        """
        result = dataframe.loc[dataframe["diffMinutes"]
                               < max_diff_minutes_for_durations].copy()
        return result

    return remove_data_under_certain_hours


def within_diff_quantile(duration_processes_quantile_threshold):
    """
    Args:
        duration_processes_quantile_threshold (float): the quantiles/percentage
        of process duration legth to keep.
    """
    def remove_data_under_quantile_for_each_event(dataframe):
        """
        Args:
            dataframe (pd.DataFrame): processed events dataframe.

        Returns:
            pd.DataFrame: processed events dataframe with extreme timings
            removed.
        """

        end_data_frame = pd.DataFrame()
        for event in dataframe["EventName"].unique():
            event_data = dataframe.loc[dataframe["EventName"] == event].copy()
            quantile_threshold = (event_data["diffMinutes"]
                               .quantile(duration_processes_quantile_threshold))
            result = event_data.loc[event_data["diffMinutes"]
                                    < quantile_threshold].copy()
            end_data_frame = pd.concat([end_data_frame, result])

        return end_data_frame

    return remove_data_under_quantile_for_each_event


def exclude_unknown_staff(dataframe):
    """
    Args:
        dataframe (pd.DataFrame): dataframe with EventStaffId column.

    Returns:
        pd.DataFrame: dataframe with staff Id 1 removed.
    """
    result = dataframe.loc[dataframe["EventStaffId"] != 1].copy()
    return result


def exclude_patients_with_uncommon_transitions_below_threshold(threshold):
    """
    Args:
        threshold (float): percentage threshold to filter out transitions.
    """
    def exclude_uncommon_transitions(dataframe):
        """
        Args:
            dataframe (pd.DataFrame): clensed events dataframe.

        Returns:
            pd.DataFrame: filtered dataframe with transfers below the threshold
            removed.
        """
        indexes_to_exclude = dataframe.loc[(dataframe["EventName"] != SPAWN)
                                           & (dataframe["EventName"] != REMOVED)
                                           & (dataframe["Percentage"] < threshold)
                                           & (dataframe["Next Event (Pathway)"]
                                              .notnull())].copy()
        visit_ids_to_exclude = indexes_to_exclude["VisitId"].unique()
        result = dataframe.loc[~dataframe["VisitId"]
                               .isin(visit_ids_to_exclude)].copy()
        return result

    return exclude_uncommon_transitions
