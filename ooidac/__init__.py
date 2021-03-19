
import numpy as np
from scipy import signal

with open("version") as fid:
    __version__ = fid.read()


def clean_dataset(dataset):
    """Remove any row in dataset for which one or more columns is np.nan
    """

    # Get rid of NaNs
    dataset = dataset[~np.isnan(dataset[:, 1:]).any(axis=1), :]

    return dataset


def boxcar_smooth_dataset(dataset, window_size):
    window = signal.windows.boxcar(window_size)
    return signal.convolve(dataset, window, 'same') / window_size
