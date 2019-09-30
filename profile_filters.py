"""

"""
import numpy as np
import logging
import os
from ooidac import processing
from configuration import DATA_CONFIG_LIST, TIMESENSOR
logger = logging.getLogger(os.path.basename(__name__))


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
        any_data = np.all(np.isnan(data))
        # if there isn't any CTD pressure data at all, we don't want the profile
        if scidata_sensor == 'sci_water_pressure' and any_data:
            remove_profile = True
            break
        allbad_scidata.append(any_data)

    if not remove_profile:
        remove_profile = np.all(allbad_scidata)

    return remove_profile


def filter_small_data_ratio(profile_data, threshold=.1, data_pts_threshold=4):
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
    timestamps = profile_data.getdata(TIMESENSOR)
    total_profile_time = timestamps[-1] - timestamps[0]

    for scidata_sensor in DATA_CONFIG_LIST:
        data = profile_data.getdata(scidata_sensor)

        # data_ratio uses ratio of data record time vs total profile time
        finites = np.flatnonzero(np.isfinite(data))
        if len(finites) >= data_pts_threshold:
            good_data_length = cum_data_time_sum(timestamps[finites])
        else:
            good_data_length = 0

        data_ratio = good_data_length/total_profile_time
        data_ratios_too_small.append(data_ratio < threshold)

    if np.all(data_ratios_too_small):
        remove_profile = True

    return remove_profile


def filter_time_lessthan(profile_data, threshold=1):
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

    minutes_of_profile = (time2 - time1) / 60.

    if minutes_of_profile < threshold:
        remove_profile = True

    return remove_profile


def filter_datatime_lessthan(profile_data, threshold=1, data_pts_threshold=4):
    """Profile filter that will remove a profile if the elapsed time for
    the data collected in a profile is less than `threshold` minutes.

    Note: a profile not removed by this filter might still be removed by
    another active filter.
    :return:
    """
    remove_profile = False
    timestamps = profile_data.getdata(TIMESENSOR)
    data_indices = processing.all_sci_indices(profile_data)

    if len(data_indices) < data_pts_threshold:
        remove_profile = True

    sci_time = timestamps[data_indices]

    minutes_of_data = cum_data_time_sum(sci_time) / 60.

    if minutes_of_data < threshold:
        remove_profile = True

    return remove_profile


def filter_no_data_at_profile_start(profile_data):
    """ Profile filter that will remove a profile if there is no science data at
     the beginning (defined as the first 10%) of the profile with extra
     emphasis on pressure from the CTD.  This tries to eliminate the case where
     a short profile that should have no data, turns on data sampling just
     before the inflection point and has enough science data to pass the other
     filters.  The rationale being that if there is no data at the start of a
     profile, it was not intended to be sampled.

    Note: a profile not removed by this filter might still be removed by
    another active filter.
    :param profile_data:
    :return: bool value if the profile is to be removed or not.
    """
    remove_profile = False
    timestamps = profile_data.getdata(TIMESENSOR)
    # ToDo: change explicit pressure here to a PRESSURESENSOR variable
    pres = profile_data.getdata('sci_water_pressure')
    first_portion_of_dive = list(range(int(len(timestamps)/10)))
    data_indices = processing.all_sci_indices(profile_data)
    pressure_ii = np.flatnonzero(np.isfinite(pres))
    if len(np.intersect1d(pressure_ii, first_portion_of_dive)) == 0:
        remove_profile = True
    if len(np.intersect1d(data_indices, first_portion_of_dive)) == 0:
        remove_profile = True

    return remove_profile


def cum_data_time_sum(sci_timestamps):
    """To eliminate the case where a small amount of science data points are at
    the beginning of a profile, and a small amount exists at the end of a
    profile, with a large non-data gap in between, this function calculates
    the cumulative sum of data time excluding time gaps larger than 3 median
    time steps.

    :param sci_timestamps: timestamps of non-nan science data records
    :return: The cumulative sum of data sample time excluding large gaps
    """
    sci_dt = np.diff(sci_timestamps)
    sci_dt_median = np.nanmedian(sci_dt)
    no_gaps_ii = sci_dt < 3 * sci_dt_median
    cum_sci_sample_time = np.sum(sci_dt[no_gaps_ii])
    return cum_sci_sample_time
