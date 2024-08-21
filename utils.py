from pathlib import Path
import pandas as pd

path_to_read_data = Path(r"./Events Data")

def sort_events(events):
    """
    Args:
        events (pd.DataFrame): clensed events dataframe.

    Returns:
        pd.DataFrame: sorted clensed events dataframe.
    """
    events = events.sort_values(by=["VisitId", "EventTime", "natural_order"])
    return events


def load_data(filename):
    """
    Args:
        filename (str): filename of the data to read in.

    Returns:
        pd.DataFrame: a dataframe of that read in data.
    """
    return pd.read_csv(f"{path_to_read_data}/{filename}")
