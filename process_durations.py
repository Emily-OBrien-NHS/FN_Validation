import math
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
import numpy as np

################################################################################
#-------------------supporting process durations functions---------------------#
################################################################################

def add_difference_in_minutes_to_durations(events_quality,
                                           where_duration_should_be_0):
    #Function to add in the time difference between each event
    #----
    # Take a copy before processing and exclude unknown staff.
    event_diffs = events_quality.copy()
    event_diffs = event_diffs.loc[event_diffs["EventStaffId"] != 1].copy()
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

def within_threshold_diff(max_diff_minutes_for_durations):
    def remove_data_under_certain_hours(dataframe):
        #Double function to allow filtering of tasks that have long time due to
        #staff member going off shift/on break if this filtering is required.
        result = dataframe.loc[dataframe["diffMinutes"]
                               < max_diff_minutes_for_durations].copy()
        return result
    return remove_data_under_certain_hours

def within_diff_quantile(duration_processes_quantile_threshold):
    def remove_data_under_quantile_for_each_event(dataframe):
        #Double function to remove extreme timings based on a specified quantile
        end_data = []
        for event in dataframe["EventName"].unique():
            event_data = dataframe.loc[dataframe["EventName"] == event].copy()
            quantile_threshold = (event_data["diffMinutes"]
                                  .quantile(
                                      duration_processes_quantile_threshold))
            result = event_data.loc[event_data["diffMinutes"]
                                    < quantile_threshold].copy()
            end_data.append(result)
        end_data_frame = pd.concat(end_data)
        return end_data_frame
    return remove_data_under_quantile_for_each_event


def only_daytime_events(dataframe):
    #Function to filter event timings to only those during 'daytime'
    eight_am_to_10_pm = dataframe.copy()
    eight_am_to_10_pm.index = dataframe["EventTime"]
    eight_am_to_10_pm = eight_am_to_10_pm.between_time("08:00:00", "20:00:00")
    result = dataframe.loc[dataframe["EventTime"]
                           .isin(eight_am_to_10_pm.index)].copy()
    return result

def generate_and_output_process_durations_log_normal(directory_path,
                                                     processed_events,
                                                     processes):
    #Function to fit a log normal distribution to the data
    #----
    #fit a log normal to each process data
    new_entries = []
    for process in processes:
        #get the data for that process
        data = (processed_events
                .loc[processed_events["Event (Pathway)"] == str(process),
                     "diffMinutes"].copy().dropna().astype(float))
        if process != "" and len(data) > 0:
            #if data for that process, fit a log normal and record parameters
            shape, loc, scale = sp.stats.lognorm.fit(data.values, floc=-0.00001)
            mu = np.log(scale)
            sigma = shape
            mean = math.exp(mu + (0.5 * sigma**2))
            variance = (math.exp(sigma**2)-1) * math.exp((2*mu)+(sigma**2))
            new_entries.append({"Event (Pathway)": str(process),
                                "Mean": mean,
                                "StdDev": math.sqrt(variance),
                                "Min": min(data.values),
                                "Max": max(data.values)})
        elif process != "":
            #if no data, record empty parameters
            new_entries.append({"Event (Pathway)": process,
                                "Mean": 0,
                                "StdDev": None,
                                "Min": None,
                                "Max": None})
    #Save lognormal results as dataframe, tidy and export to csv
    process_durations = pd.DataFrame(new_entries)
    process_durations["Notes"] = None
    process_durations = process_durations.rename(
                        columns={"Mean": "Duration Mean",
                        "Event (Pathway)": "Process (Pathway and Recurrent)"})
    process_durations.to_csv(directory_path/"Process Durations.csv",index=False)

################################################################################
#----------------main histogram and process durations function-----------------#
################################################################################

def main_generate_histogram_and_process_durations(directory_path,
    processed_events, groupby_column, output_path, plots, filterFuncs=None):
    #Main function to fit lognormals, plot distributions if required and output
    #the process durations input file.
    #----
    #Make output folder.
    plot_folder_directory_path = output_path / "Durations" / directory_path
    plot_folder_directory_path.mkdir(exist_ok=True, parents=True)
    #Appply functions if these have been passed in.
    if filterFuncs is not None:
        for filter in filterFuncs:
            processed_events = filter(processed_events)
    #Get the average number of Treatment repeats and multiply the treatment
    #time by this to create 'mega treatment' event.
    treat_repeat = ((processed_events
                     .loc[processed_events['EventName'] == 'Treatment']
                     .groupby(['Pathway', 'VisitId'], as_index=False)
                      ['EventTime'].count()
                      .groupby('Pathway')['EventTime'].mean())
                     .rename('TreatRepeat'))
    processed_events = processed_events.merge(treat_repeat, on='Pathway')
    processed_events['diffMinutes'] = (
        np.where(processed_events['EventName'] == 'Treatment',
                 processed_events['diffMinutes']
                 * processed_events['TreatRepeat'],
                 processed_events['diffMinutes']))
    processed_events = processed_events.drop('TreatRepeat', axis=1)
    #Log normal distributions for each process
    processes = processed_events["Event (Pathway)"].unique().tolist()
    generate_and_output_process_durations_log_normal(plot_folder_directory_path,
                                                    processed_events, processes)
    if plots:
        event_times = (processed_events.groupby(groupby_column)["diffMinutes"])
        for group, data in event_times:
            #create a histogram for the data in each event.
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
