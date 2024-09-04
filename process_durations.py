
import math
import pandas as pd
import scipy as sp  # type: ignore
import matplotlib.pyplot as plt
import numpy as np
from event_filtering_functions import exclude_unknown_staff
from config import where_duration_should_be_0


def generate_and_output_process_durations_log_normal(directory_path,
                                                    processed_events, processes):
    """
    Args:
        directory_path (Path): path to output directory.
        processed_events (pd.DataFrame): dataframe of processed events.
        processes (list[str]): list of processes.
    """
    #fit a log normal to each process data
    new_entries = []
    for process in processes:
        #get the data for that process
        data = (processed_events.loc[processed_events["Event (Pathway)"]
                                    == str(process), "diffMinutes"].copy().dropna()
                                    .astype(float))
        if process != "" and len(data) > 0:
            #if data for that process, fit a log normal and record parameters
            shape, loc, scale = sp.stats.lognorm.fit(data.values, floc=-0.00001)
            mu = np.log(scale)
            sigma = shape
            mean = math.exp(mu + (0.5 * sigma**2))
            variance = (math.exp(sigma**2) - 1) * math.exp((2 * mu) + (sigma**2))
            new_entries.append({"Event (Pathway)": str(process),
                                          "Mean": mean,
                                          "StdDev": math.sqrt(variance),
                                          "Min": min(data.values),
                                          "Max": max(data.values)})
        elif process != "":
            #if no data, record empty parameters
            new_entries.append({"Event (Pathway)": process, "Mean": 0,
                                "StdDev": None, "Min": None, "Max": None})

    #Save lognormal results as dataframe, tidy and export to csv
    process_durations = pd.DataFrame(new_entries)
    process_durations["Notes"] = None
    process_durations = process_durations.rename(columns={
        "Mean": "Duration Mean",
        "Event (Pathway)": "Process (Pathway and Recurrent)"})
    process_durations.to_csv(directory_path / "Process Durations.csv",
                             index=False)


def generate_and_output_histogram_and_process_durations(directory_path,
    processed_events, groupby_column, output_path, plots, filterFuncs=None):
    """
    Args:
        directory_path (str): path to create/navigate to output directory.
        processed_events (pd.DataFrame): dataframe of processed events.
        groupby_column (str): column name to group by.
        output_path (Path): path to output folder.
        filterFuncs (_type_, optional): Functions to apply to data if required.
        Defaults to None.
    """
    #Make output folder.
    plot_folder_directory_path = output_path / "Durations" / directory_path
    plot_folder_directory_path.mkdir(exist_ok=True, parents=True)

    #Get a list of unique events.
    processes = processed_events["Event (Pathway)"].unique().tolist()

    #Appply functions if these have been passed in.
    if filterFuncs is not None:
        for filter in filterFuncs:
            processed_events = filter(processed_events)

    #Get the average number of Treatment repeats and multiply the treatment
    #time by this to create 'mega treatment' event.
    treat_repeat = (processed_events.loc[processed_events['EventName'] == 'Treatment']
                    .groupby(['Pathway', 'VisitId'], as_index=False)['EventTime']
                    .count().groupby('Pathway')['EventTime'].mean()).rename('TreatRepeat')
    processed_events = processed_events.merge(treat_repeat, on='Pathway')
    processed_events['diffMinutes'] = np.where(
                                      processed_events['EventName'] == 'Treatment',
                                      processed_events['diffMinutes']
                                       * processed_events['TreatRepeat'],
                                      processed_events['diffMinutes'])
    processed_events = processed_events.drop('TreatRepeat', axis=1)

    #Log normal distributions.
    generate_and_output_process_durations_log_normal(plot_folder_directory_path,
                                                    processed_events, processes)
    if plots:
        for group, data in processed_events.groupby(groupby_column)["diffMinutes"]:
            #create a histogram for the data in each event.
            if 'Treatment' in group:
                data = data *1.3
            title = str(group).replace("_", "")
            ax = data.plot.hist(bins=100, figsize=(12, 8))
            ax.grid(True, which="both", linestyle="--", linewidth=0.5)
            ax.set_xlabel("time (minutes)")
            ax.set_title(f"{title} {directory_path}")
            fig = ax.get_figure()
            if fig is not None:
                slash_replacement = {"\\": "-", "//": "-", "/": "-"}
                for key, value in slash_replacement.items():
                    label = f"{(str(title)).replace(key, value)}.png"
                fig.savefig(plot_folder_directory_path / f"{label}")
                plt.close(fig)


def add_difference_in_minutes_to_durations(events_quality):
    """
    Args:
        events_quality (pd.DataFrame): clensed events dataframe.

    Returns:
        pd.DataFrame: events dataframe with durations added.
    """
    # Take a copy before processing and exclude unknown staff.
    event_diffs = events_quality.copy()
    event_diffs = exclude_unknown_staff(event_diffs)

    # Re-sort to ensure in DateTime order
    event_diffs = event_diffs.sort_values(by=["EventStaffId", "EventTime"])

    # get time diff column
    event_diffs["time_diff"] = (event_diffs.groupby("EventStaffId")
                                ["EventTime"].diff().shift(-1))
    event_diffs["diffMinutes"] = (pd.to_timedelta(event_diffs["time_diff"])
                                  .dt.total_seconds() / 60)

    event_diffs.loc[event_diffs["EventName"].isin(where_duration_should_be_0)
                    & event_diffs["diffMinutes"].isna(), "diffMinutes"] = 0
    return event_diffs
