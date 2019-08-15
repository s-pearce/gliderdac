import numpy as np
import matplotlib.pyplot as plt
import datetime as dt


def fwd_fill(x, init_nan_fill=None):
    """ Forward fills an array with nans, i.e. fill each nan or nan cluster with
    the preceding finite value.

    :param x: array with nan clusters to forward fill
    :param init_nan_fill: value to fill with if there are initial nans,
        defaults to back filling with the first finite value.
    :return: array x with nans forward filled
    """

    # This function assumes that the first value in the array is finite.
    # That should be accounted for.  This can be accounted for by identifying
    # nan category number 0. if there is a zero, it means it is at the start
    # of the array.
    y = x.copy()

    # first get indexes of nans
    nan_index = np.where(np.isnan(y))[0]  # assuming 1d, so take [0] index

    # create an ordinal set vector the same length as nan_index
    # i.e. (0,1,2,3,4,...)
    nanrange = range(len(nan_index))

    # subtraction basically makes a numerical category for each consecutive
    # index value

    # subtraction creates a "category" vector where each unique value
    # represents a specific sequence of nans in y, and, the number of
    # occurrences of a particular value gives the number of NaNs in that
    # particular sequence.

    # E.g.
    # x = [19.8, nan, nan, nan, 11.3, 9.3, nan, nan, nan, 7.4, nan, nan, 4.5])
    # nan_index = [1, 2, 3, 6, 7, 8, 10, 11]
    # nanrange =  [0, 1, 2, 3, 4, 5,  6,  7]; subtract nanrange from nan_index
    # nan_cats =  [1, 1, 1, 3, 3, 3,  4,  4];
    # categorical numbers represent a group of consecutive nans
    nan_cats = nan_index - nanrange

    # using np.unique creates the unique category groups, and also returns an
    # index where the first occurrence of the category group number occurs (
    # firstidx), and an array of indexes that can reproduce nan_cats from
    # uniques (restoreidx)
    uniques, firstidx, restoreidx = np.unique(nan_cats, True, True)

    # using firstidx on nan_index instead, gives the indexes of the first nan
    # in each cluster.  Subtracting those indexs by 1 gives the indexes for
    # the finite value preceeding each nan cluster.  Then simply use the
    # restoring index to repeat the preceeding finite values into an index
    # array (good_index) that is the same shape as nan_index.
    good_index = (nan_index[firstidx] - 1)[restoreidx]

    # So using x[good_index] makes an array of the preceeding good values with
    # the same repititions as the nan clusters.  Then you can reinsert that
    # repeating finite value array into where the nans were in the original
    # array using ...

    if not init_nan_fill:
        init_nan_fill = y[np.isfinite(y)][0]
    if 0 in uniques:
        y = np.append(y, init_nan_fill)

    y[nan_index] = y[good_index]

    if 0 in uniques:
        y = np.delete(y, -1)

    return y


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
