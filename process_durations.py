import math
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
import numpy as np
import config

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

def add_imaging_timings(event_diffs, imaging_raw):
    #Remove imaging data from the events file
    event_diffs = event_diffs.loc[~event_diffs["EventName"]
                                  .isin(config.imaging_events)].copy()
    #Filter imaging to MRI and Ultrasound, and get their timings
    image_times = imaging_raw.loc[imaging_raw["EventName"]
                                  .isin(["MRI", "Ultrasound"])].copy()
    image_times["diffMinutes"] = (pd.to_timedelta(
        image_times["ResultsAvailableDateTime"] - image_times["EventTime"])
        .dt.total_seconds() / 60)
    image_times["Event (Pathway)"] = image_times["EventName"]
    #Concat this back onto the events file
    event_diffs = pd.concat([event_diffs,image_times[["VisitId", "EventTime",
                            "EventName", "Event (Pathway)", "diffMinutes"]]])
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

def plot_distribution(process, data, shape, scale, plot_path):
    #If plots is true, create the plot of the histogram and distribution
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    #got log normal and an array of x values to plot it
    Y = sp.stats.lognorm(s=shape, scale=scale)
    x = np.linspace(min(data), max(data), 1000)
    #plot histogram
    ax.hist(data, bins=100, density=True, label=process)
    #plt log normal
    ax.plot(x, Y.pdf(x), label='lognormal', color='red')
    #add trimmings and save
    ax.legend()
    ax.set_xlabel("time (minutes)")
    file_name = process.replace('/', '-')
    ax.set_title(process)
    fig.savefig(plot_path + f"/{file_name}")
    plt.close()

def timings_dict(process, mean, std, min_, max_, note):
    new_entries_dict = {"Event (Pathway)":process,
                       "Mean": mean,
                       "StdDev": std,
                       "Min": min_,
                       "Max": max_,
                       "Notes":note}
    return new_entries_dict

def generate_and_output_process_durations_log_normal(directory_path, plot_path,
                                                     processed_events,
                                                     processes):
    #Function to fit a log normal distribution to the data
    #----
    #fit a log normal to each process data
    new_entries = []
    for process in processes:
        #get the data for that process, with 0s removed
        data = (processed_events.loc[
                           (processed_events["Event (Pathway)"] == str(process))
                           & (processed_events['diffMinutes'] > 0),
                           "diffMinutes"].copy().dropna().astype(float))
        if process != "" and len(data) > 0:
            #if data for that process, fit a log normal (use mean of nlog of
            # data as scale start point) and record parameters
            shape, loc, scale = sp.stats.lognorm.fit(
                          data.values, scale=np.log(data).mean(), floc=-0.00001)
            mu = np.log(scale)
            sigma = shape
            mean = math.exp(mu + (0.5 * sigma**2))
            variance = (math.exp(sigma**2)-1) * math.exp((2*mu)+(sigma**2))
            std = math.sqrt(variance)
            min_ = min(data.values)
            max_ = max(data.values)
            #If process not on a pathway, ensure this gets added later.
            pathway_loop = not "(" in process
            #plot the distribution
            if config.plots:
                plot_distribution(process, data, shape, scale, plot_path)
        else:
            #If no data, add 0s and np.nan
            mean = 0
            std, min_, max_ = np.nan
            pathway_loop = False

        #If process is not on a pathway, repeat the timings for each pathway,
        #otherwise just add the data.
        if pathway_loop:
            for pathway in config.pathways:
                new_entries.append(
                    timings_dict((str(process) + " (" + pathway + ")"),
                                 mean, std, min_, max_, np.nan))
        else:
            new_entries.append(timings_dict(process, mean, std, min_, max_,
                                            np.nan))

    #----Add in manually added process timings
    for process, times in config.add_process_durs.items():
        mean, std, min_, max_ = times
        #If event doesn't belong to a pathway, add this in
        if '(' not in process:
            for pathway in config.pathways:
                new_entries.append(
                    timings_dict((str(process) + " (" + pathway + ")"),
                                 mean, std, min_, max_, "Timings from config"))
        else:
            new_entries.append(timings_dict(process, mean, std, min_, max_,
                                            "Timings from config"))

    #----Save lognormal results as dataframe, tidy and export to csv
    process_durations = pd.DataFrame(new_entries)
    process_durations = process_durations.rename(
                        columns={"Mean": "Duration Mean",
                        "Event (Pathway)": "Process (Pathway and Recurrent)"})
    process_durations.to_csv(directory_path + "/Process Durations.csv",index=False)

################################################################################
#----------------main histogram and process durations function-----------------#
################################################################################

def main_generate_histogram_and_process_durations(
    processed_events, treat_repeat, output_path, filterFuncs=None):
    #Main function to fit lognormals, plot distributions if required and output
    #the process durations input file.
    #----
    #Specify output folder.
    plot_folder_path = output_path + "/Additional Outputs/Duration Distributions"
    #Appply functions if these have been passed in.
    if filterFuncs is not None:
        for filter in filterFuncs:
            processed_events = filter(processed_events)
    #multiply the treatment time by the average number of treatment events per
    #patient to create a 'mega treatment' event.
    processed_events = processed_events.merge(treat_repeat, on='Pathway', how='left')
    processed_events['diffMinutes'] = (
        np.where(processed_events['EventName'] == 'Treatment',
                 processed_events['diffMinutes']
                 * processed_events['TreatRepeat'],
                 processed_events['diffMinutes']))
    processed_events = processed_events.drop('TreatRepeat', axis=1)
    #Log normal distributions for each process
    processes = processed_events["Event (Pathway)"].unique().tolist()
    generate_and_output_process_durations_log_normal(output_path, plot_folder_path,
                                                     processed_events, processes)
