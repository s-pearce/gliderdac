"""

"""
import numpy as np
import logging
import os
from configuration import DATA_CONFIG_LIST
logger = logging.getLogger(os.path.basename(__name__))

TIMESENSOR = 'm_present_time'


def filter_no_data(profile_data):
    """Profile filter that will remove a profile if all of the relevant science
    sensors' data (listed by the SCI_DATA_PROFILE_LIST configuration parameter)
    is missing.  It keeps the profile if ANY of the science sensors have any
    data present.

    Note: a profile not removed by this filter might still be removed by
    another active filter.
    :return:
    """
    remove_profile = False
    allbad_scidata = []

    for scidata_sensor in DATA_CONFIG_LIST:
        data = profile_data.getdata(scidata_sensor)
        allbad_scidata.append(np.all(np.isnan(data)))

    if np.all(allbad_scidata):
        remove_profile = True

    return remove_profile


def filter_small_data_ratio(profile_data, threshold=.01):
    """Profile filter that will remove a profile if all of the relevant science
    sensors ( listed by the SCI_DATA_PROFILE_LIST configuration parameter)
    have a ratio of good data to missing data that is smaller than the
    threshold.  It keeps the profile if ANY of the science sensors have a
    ratio larger than the threshold.

    Note: a profile not removed by this filter might still be removed by
    another active filter.
    :return:
    """
    remove_profile = False
    data_ratios_too_small = []
    for scidata_sensor in DATA_CONFIG_LIST:
        data = profile_data.getdata(scidata_sensor)
        timestamps = profile_data.getdata(TIMESENSOR)

        # data_ratio version 1 uses points
        # good_data_length = len(np.flatnonzero(np.isfinite(data)))
        # total_data_length = len(data)

        # data_ratio version 2 uses time length
        total_data_length = timestamps[-1] - timestamps[0]
        finites = np.flatnonzero(np.isfinite(data))
        if len(finites) > 0:
            good_data_length = timestamps[finites[-1]] - timestamps[finites[0]]
        else:
            good_data_length = 0

        data_ratio = good_data_length/total_data_length
        data_ratios_too_small.append(data_ratio < threshold)

    if np.all(data_ratios_too_small):
        remove_profile = True

    return remove_profile


def filter_time_lessthan(profile_data, threshold=5):
    """Profile filter that will remove a profile if the elapsed time for
    the profile is less than `threshold` minutes.

    Note: a profile not removed by this filter might still be removed by
    another active filter.
    :return:
    """
    remove_profile = False
    timestamps = profile_data.getdata(TIMESENSOR)
    time1 = timestamps[0]
    time2 = timestamps[-1]

    minutes_of_data = (time2 - time1) / 60.

    if minutes_of_data < threshold:
        remove_profile = True

    return remove_profile
