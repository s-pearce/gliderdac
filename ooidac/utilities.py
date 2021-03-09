import numpy as np
import matplotlib.pyplot as plt
import datetime as dt


def cluster_index(indices, ids=False):
    ordinals = np.arange(0, len(indices))
    index = indices - ordinals
    if ids:
        index_ids = np.unique(index)
        return index, index_ids
    return index


def nan_array_equal(arr1, arr2):
    finites1 = np.isfinite(arr1)
    finites2 = np.isfinite(arr2)
    aa = np.array_equal(finites1, finites2)

    nans1 = np.isnan(arr1)
    nans2 = np.isnan(arr2)
    bb = np.array_equal(nans1, nans2)

    cc = np.array_equal(arr1[finites1], arr2[finites2])

    if aa and bb and cc:
        return True
    else:
        return False


def plot_vert_lines(x, style='k:', ax=None):
    x = np.atleast_1d(x)
    n = len(x)
    x = np.tile(x, (2, 1))
    if not ax:
        ax = plt.gca()
    ylims = ax.get_ylim()
    y = np.array([np.repeat(ylims[0], n), np.repeat(ylims[1], n)])
    plt.plot(x, y, style)


def date_xlabel(time_range):
    """date_xlabel
    For use with plotting time on the x-axis, create_date_label provides a
    overview date range for the xlabel of the plot.  Current implementation
    is when the plotted data is over the range of hours to days and assumes
    that matplotlib auto generated meaningful time xtick marks from a
    datetime object.

    :param time_range: sequence of datetime objects in sequential order.
    """
    if len(time_range) < 2:
        raise InputError("time_range must be a sequence of date times")
    t1 = time_range[0]
    t2 = time_range[-1]
    if not isinstance(t1, dt.datetime) or not isinstance(t2, dt.datetime):
        raise InputError(
            "Time Range sequence elements must be datetime objects")

    date_label = t1.date().isoformat()

    year_diff = t1.year != t2.year
    mnth_diff = t1.month != t2.month
    day_diff = t1.day != t2.day
    if year_diff or mnth_diff or day_diff:
        date_label += (
            "--" +
            (year_diff * ("{:02d}".format(t2.year) + "-")) +
            (mnth_diff * ("{:02d}".format(t2.month) + "-")) +
            (day_diff * "{:02d}".format(t2.day))
        )
    return date_label


class InputError(Exception):
    pass
