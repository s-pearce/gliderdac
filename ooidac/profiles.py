#!/usr/bin/env python

import numpy as np
import logging
import os

from ooidac import boxcar_smooth_dataset

from ooidac.utilities import fwd_fill
from ooidac.processing import all_sci_indices
import profile_filters
# import pdb

# ToDo: fix the imports above

logger = logging.getLogger(os.path.basename(__file__))


class Profiles(object):
    def __init__(self, dba):
        """ Given a glider DbaData object, a class for storing and
        spliting the glider dives into profiles (I.e. a sampling dive or
        climb).

        :param dba: DbaData object of ascii glider data with the llat
            variables added
        """
        self.dba = dba
        self._indices = None
        self.inflection_times = None
        # Todo: add in some checking before proceeding.

    def __iter__(self):
        """makes it possible to iterate over the profiles in a for loop by
        returning a GliderData object that is sliced by the profile indices"""
        if self._indices is not None:
            for profile_indices in self._indices:
                yield self.dba.slicedata(indices=profile_indices)

    def __len__(self):
        """allows the len function to be called on a Profiles instance and it
        will return the number of profiles found"""
        if self.indices:
            return len(self.indices)
        else:
            return 0

    @property
    def indices(self):
        return self._indices

    def get_profile(self, index):
        try:
            return self.dba.slicedata(indices=self._indices[index])
        except IndexError:
            raise ProfileIndexError("Profile index is out of range")

    # TODO: fix all of the return statements in this class (maybe?), I think
    #  I meant here that sometimes the return statements return an empty
    #  list, when they could just say "return" alone.

    def find_profiles_by_depth(
            self, depth_sensor='m_depth', tsint=2, winsize=10):
        """Discovery of profiles in a glider segment using depth and time.

        Profiles are discovered by smoothing the depth timeseries and using the
        derivative of depth vs time to find the inflection points to break
        the segment into profiles.  Profiles are truncated to where science
        data exists; science sensors for the glider are configured in
        Configuration.py.
        Depth can be `m_depth` or CTD pressure (Default m_depth).  For smoothing
        a filtered depth is created at regular time intervals `tsint` (default
        2 secs) and boxcar filtered with window `winsize` (default 5 points).
        The smoothing is affected by both choice of `tsint` and `winsize`, but
        usually still returns similar profiles.  After profiles are discovered,
        they should be filtered with the `filter_profiles` method of this
        class, which removes profiles that are not true profiles.

        :param depth_sensor: The Depth sensor to use for profile discovery.
        Should be either m_depth, sci_water_pressure, or a derivative of
        those 2.  Default is `m_depth`
        :param tsint: Time interval in seconds for filtered depth.
        This affects filtering.  Default is 2.
        :param winsize: Window size for boxcar smoothing filter.
        :return: output is a list of profile indices in self.indices
        """
        self._indices = []
        depth = self.dba.getdata(depth_sensor)
        time_ = self.dba.getdata('m_present_time')

        # Set negative depth values to NaN; using this method to avoid numpy
        # warnings from < when nans are in the array
        depth_ii = np.flatnonzero(np.isfinite(depth))  # non-nan indices
        neg_depths = np.flatnonzero(depth[depth_ii] <= 0)  # indices to depth_ii
        depth[depth_ii[neg_depths]] = np.nan

        # Remove NaN depths and truncate to when science data begins being
        # recorded and ends
        depth_ii = np.flatnonzero(np.isfinite(depth))

        sci_indices = all_sci_indices(self.dba)  # from ooidac.processing
        if len(sci_indices) > 0:
            starting_index = sci_indices[0]
            ending_index = sci_indices[-1]
        else:
            return  # no science_indices, then we don't care to finish

        depth_ii = depth_ii[
            np.logical_and(
                depth_ii >= starting_index,
                depth_ii <= ending_index)
        ]

        # ---Create a smoothed depth timeseries for finding inflections ------#

        # Find start and end times first adding winsize * tsint timesteps
        # onto the start and end to account for filter edge effects
        itime_start = np.ceil(time_[depth_ii].min()) - winsize * tsint
        itime_end = np.floor(time_[depth_ii].max()) + (winsize + 1) * tsint

        itime = np.arange(itime_start, itime_end, tsint)
        idepth = np.interp(itime, time_[depth_ii], depth[depth_ii],
                           left=depth[depth_ii[0]], right=depth[depth_ii[-1]])
        fz = boxcar_smooth_dataset(idepth, winsize)

        # remove the extra points with filter edge effects
        fz = fz[winsize:-winsize]
        itime = itime[winsize:-winsize]
        idepth = idepth[winsize:-winsize]

        # Zero crossings of the time derivative of filtered depth are the
        # inflection points.  Differential time is midway between the
        # filtered timestamps.
        # Originally, this used scipy's fsolver to locate the exact zero
        # crossing, but only the timestamp before the zero crossing is needed
        # to be the last in a profile and the timestamp after the zero
        # crossing to be the first in the next profile.

        dz_dt = np.diff(fz) / np.diff(itime)
        dtime = itime[:-1] + np.diff(itime) / 2  # differential time

        # Get the time point just after a zero crossing.  The flatnonzero
        # statement below gets the point before a zero crossing.
        zero_crossings_ii = np.flatnonzero(abs(np.diff(np.sign(dz_dt))))
        zc_times = dtime[zero_crossings_ii] + (
                dtime[zero_crossings_ii + 1] - dtime[zero_crossings_ii]) / 2.

        profile_switch_times = zc_times[np.logical_and(
            zc_times > time_[starting_index],
            zc_times < time_[ending_index]
        )]
        # insert the timestamp of the first science data point at the start
        # and the last data point at the end.
        profile_switch_times = np.insert(
            profile_switch_times, [0, len(profile_switch_times)],
            [time_[starting_index],
             time_[ending_index]])

        self.inflection_times = profile_switch_times

        # profile_switch_times = self.adjust_inflections(depth, time_)
        profile_switch_times = self.adjust_inflections(depth, time_)

        # use the time range to gather indices for each profile
        for ii in range(len(profile_switch_times)-1):
            pstart = profile_switch_times[ii]
            pend = profile_switch_times[ii+1]
            profile_ii = np.flatnonzero(
                np.logical_and(
                    time_ >= pstart,
                    time_ <= pend))  # inclusive since before the inflection
            if len(profile_ii) == 0:
                continue
            self._indices.append(profile_ii)

    def adjust_inflections(self, depth, time_):
        """Filters out bad inflection points.

        Bad inflection points are small surface, bottom of dive, or mid-profile
        wiggles that are not associated with true dive or climb inflections.
        These false inflections are removed so that when profile indices are
        created, they don't separate into separate small profiles.

        :param depth:
        :param time_:
        :return:
        """
        inflections = self.inflection_times
        inflection_depths = np.interp(
            inflections, time_[np.isfinite(depth)],
            depth[np.isfinite(depth)]
        )

        # First remove the false diving inflections (i.e. the small wiggles) by
        # taking the good inflection and looking ahead until an inflection depth
        # difference greater than 2m is found
        inflx_ii = 0
        fwd_counter = 1
        inflx_to_keep = np.full(len(inflections), True)
        while inflx_ii < len(inflections):
            ii_depth = inflection_depths[inflx_ii]  # depth of current inflection
            # look ahead for the next true inflection change
            if inflx_ii + fwd_counter >= len(inflections):
                break
            while abs(inflection_depths[inflx_ii + fwd_counter] - ii_depth) < 2:
                inflx_to_keep[inflx_ii + fwd_counter] = False
                fwd_counter += 1
                if inflx_ii + fwd_counter >= len(inflections):
                    break
            inflx_ii = inflx_ii + fwd_counter
            fwd_counter = 1

        # afterwards we may be left with mid profile direction changes that were
        # greater than 2 m.  But now they can identified by not changing trend,
        # since any of our good  inflection points left will change trend sign.
        good_inflx_ii = np.flatnonzero(inflx_to_keep)
        trends = np.diff(inflection_depths[good_inflx_ii])
        # find where the trends are the same by getting the diff of the sign of
        # the trend (which is also a diff).  Must add one because diff always
        # results in N-1
        same_trends = np.flatnonzero(np.diff(np.sign(trends)) == 0) + 1
        good_inflx_ii = np.delete(good_inflx_ii, same_trends)
        self.inflection_times = inflections[good_inflx_ii]

        return inflections[good_inflx_ii]

    def adjust_inflections_old(self, depth, time_):
        """Filters out bad inflection points.

        Bad inflection points are small surface, bottom of dive, or mid-profile
        wiggles that are not associated with true dive or climb inflections.
        These false inflections are removed so that when profile indices are
        created, they don't separate into separate small profiles.

        :param depth:
        :param time_:
        :return:
        """
        inflection_times = self.inflection_times
        inflx_to_keep = np.full(len(inflection_times), True)
        inflection_depths = np.interp(
            inflection_times, time_[np.isfinite(depth)],
            depth[np.isfinite(depth)]
        )
        for ii in range(1, len(inflection_times) - 1):
            # for each inflection point, get the inflection point before and
            # afer to compare depths.
            zcs = inflection_times[
                  ii - 1:ii + 2]
            flx_depths = inflection_depths[ii - 1:ii + 2]

            # Taking the 2nd derivative tells us the inflection direction.
            # Positive inflection is toward larger depths (i.e.
            # toward the sea floor) and negative inflections are toward
            # smaller depths (i.e. toward the surface)
            dz_dt = np.diff(flx_depths) / np.diff(zcs)
            dt = zcs[:-1] + np.diff(zcs) / 2
            d2z_dt2 = np.diff(dz_dt) / np.diff(dt)

            # Inflections where both inflection points to either side have less
            # than 2m difference are no good
            if (
                    abs(flx_depths[1] - flx_depths[0]) < 2
                    and abs(flx_depths[1] - flx_depths[2]) < 2):
                inflx_to_keep[ii] = False
                continue

            # If a positive inflection is not the top of a yo, it is no
            # good.  This is determined by positive inflections (climb to
            # dive, i.e. near surface) where the inflections to either side
            # are less than 2m deep are no good. Positive inflections should
            # be bordered by inflections that are deep enough to consider
            # that it dove, here chosen to be anything deeper than 2 m.
            if (
                    d2z_dt2[0] > 0
                    and (flx_depths[0] < 2 or flx_depths[2] < 2)):
                inflx_to_keep[ii] = False
                continue

            # negative inflections (dive to climb, i.e. at the bottom of a
            # dive) that are less than the depth considered to be a dive (
            # here 2m) are no good.
            if d2z_dt2[0] < 0 and flx_depths[1] < 2:
                inflx_to_keep[ii] = False
                continue

            # if an inflection point's 2nd derivative is zero (meaning not
            # really inflecting) then it is no good.
            if d2z_dt2[0] == 0:
                inflx_to_keep[ii] = False

        # Then remove those inflections that are True in inflx_to_rm
        inflection_times = inflection_times[inflx_to_keep]
        inflection_depths = inflection_depths[inflx_to_keep]

        # After, we now need to see if there are any mid profile breaks/
        # false inflections (might include a small change in the wrong
        # direction) and remove those inflection points.

        # 1. To do this, start by seeing if there are any differences left
        # between inflection depths that are less less than 2m and remove
        # them (remove the one after the small difference).
        inflx_to_keep = np.full(len(inflection_times), True)
        small_diffs = np.flatnonzero(
            abs(np.diff(inflection_depths)) < 2) + 1  # +1 for after
        inflx_to_keep[small_diffs] = False
        inflection_times = inflection_times[inflx_to_keep]
        inflection_depths = inflection_depths[inflx_to_keep]

        # 2. Again take the difference, and where the trend between 2
        # sequential inflections have the same sign (i.e. they both trend in
        # the same direction) remove the second that has the same sign.
        inflx_to_keep = np.full(len(inflection_times), True)
        same_sign_ii = np.flatnonzero(
            np.diff(np.sign(np.diff(inflection_depths))) == 0)
        same_sign_ii += 1  # add one to each index since we took diff twice
        inflx_to_keep[same_sign_ii] = False
        return inflection_times[inflx_to_keep]

    def find_profiles_by_depth1(
            self, depth_sensor='m_depth', tsint=2, winsize=5):
        """Discovery of profiles in a glider segment using depth and time.

        Profiles are discovered by smoothing the depth timeseries and using the
        derivative of depth vs time to find the inflection points to break
        the segment into profiles.  Profiles are truncated to where science
        data exists; science sensors for the glider are configured in
        Configuration.py.
        Depth can be `m_depth` or CTD pressure (Default m_depth).  For smoothing
        a filtered depth is created at regular time intervals `tsint` (default
        2 secs) and boxcar filtered with window `winsize` (default 5 points).
        The smoothing is affected by both choice of `tsint` and `winsize`, but
        usually still returns similar profiles.  After profiles are discovered,
        they should be filtered with the `filter_profiles` method of this
        class, which removes profiles that are not true profiles.

        :param depth_sensor: The Depth sensor to use for profile discovery.
        Should be either m_depth, sci_water_pressure, or a derivative of
        those 2.  Default is `m_depth`
        :param tsint: Time interval in seconds for filtered depth.
        This affects filtering.  Default is 2.
        :param winsize: Window size for boxcar smoothing filter.
        :return: output is a list of profile indices in self.indices
        """
        self._indices = []
        depth = self.dba.getdata(depth_sensor)
        time_ = self.dba.getdata('m_present_time')

        # Set negative depth values to NaN; using this method to avoid numpy
        # warnings from < when nans are in the array
        depth_ii = np.flatnonzero(np.isfinite(depth))  # non-nan indices
        neg_depths = np.flatnonzero(depth[depth_ii] <= 0)  # indices to depth_ii
        depth[depth_ii[neg_depths]] = np.nan

        # Remove NaN depths and truncate to when science data begins being
        # recorded and ends
        depth_ii = np.flatnonzero(np.isfinite(depth))

        sci_indices = all_sci_indices(self.dba)  # from ooidac.processing
        if len(sci_indices) > 0:
            starting_index = sci_indices[0]
            ending_index = sci_indices[-1]
        else:
            return  # no science_indices, then we don't care to finish

        depth_ii = depth_ii[
            np.logical_and(
                depth_ii >= starting_index,
                depth_ii <= ending_index)
        ]

        # ---Create a smoothed depth timeseries for finding inflections ------#

        # find start and end times first adding winsize * tsint timesteps
        # onto the start and end to account for filter edge effects
        itime_start = np.ceil(time_[depth_ii].min()) - winsize * tsint
        itime_end = np.floor(time_[depth_ii].max()) + (winsize + 1) * tsint

        itime = np.arange(itime_start, itime_end, tsint)
        idepth = np.interp(itime, time_[depth_ii], depth[depth_ii],
                           left=depth[depth_ii[0]], right=depth[depth_ii[-1]])
        fz = boxcar_smooth_dataset(idepth, winsize)

        # remove the extra points with filter edge effects
        fz = fz[winsize:-winsize]
        itime = itime[winsize:-winsize]
        idepth = idepth[winsize:-winsize]

        # Zero crossings of the time derivative of filtered depth are the
        # inflection points.  Differential time is midway between the
        # filtered timestamps.
        # Originally, this used scipy's fsolver to locate the exact zero
        # crossing, but only the timestamp before the zero crossing is needed
        # to be the last in a profile and the timestamp after the zero
        # crossing to be the first in the next profile.

        dz_dt = np.diff(fz) / np.diff(itime)
        dtime = itime[:-1] + np.diff(itime) / 2  # differential time

        # Get the time point just after a zero crossing.  The flatnonzero
        # statement below gets the point before a zero crossing.
        zero_crossings_ii = np.flatnonzero(abs(np.diff(np.sign(dz_dt))))
        zc_times = dtime[zero_crossings_ii] + (
                dtime[zero_crossings_ii + 1] - dtime[zero_crossings_ii]) / 2.
        self.inflection_times = zc_times

        profile_switch_times = zc_times[np.logical_and(
            zc_times > time_[starting_index],
            zc_times < time_[ending_index]
        )]
        # insert the timestamp of the first science data point at the start
        # and the last data point at the end.
        profile_switch_times = np.insert(
            profile_switch_times, [0, len(profile_switch_times)],
            [time_[starting_index],
             time_[ending_index]])

        # create a time range for each profile
        profile_times = np.empty((len(profile_switch_times) - 1, 2))
        profile_times[:, 0] = profile_switch_times[:-1]
        profile_times[:, 1] = profile_switch_times[1:]

        # use the time range to gather indices for each profile
        for pstart, pend in profile_times:
            profile_ii = np.flatnonzero(
                np.logical_and(
                    time_ >= pstart,
                    time_ <= pend))  # inclusive since before the inflection
            if len(profile_ii) == 0:
                continue
            self._indices.append(profile_ii)

    def orig_find_profiles_by_depth(self, tsint=2, filter_winsize=10):
        """Returns the start and stop timestamps for every profile indexed from
        the depth timeseries

        Parameters:
            time, depth

        Returns:
            A Nx2 array of the start and stop timestamps indexed from the yo

        Use filter_yo_extrema to remove invalid/incomplete profiles
        """

        # Create list of profile indices - pearce / kerfoot method
        self._indices = []

        if 'llat_time' in self.dba.sensor_names:
            timestamps = self.dba['llat_time']
        else:
            timestamps = self.dba['m_present_time']

        timestamps = timestamps['data']
        if 'm_depth' not in self.dba.sensor_names:
            logging.warning(
                ('m_depth not found in dba {:s} for profiles, '
                 'trying pressure/depth instead').format(self.dba.source_file)
            )
            if 'llat_depth' not in self.dba.sensor_names:
                logging.warning('no depth source found in dba {:s}'.format(
                    self.dba.source_file)
                )
                return
            depth = self.dba['llat_depth']['data']
        else:
            depth = self.dba['m_depth']['data']

        # validate_glider_args(timestamps, depth)

        # Set negative depth values to NaN; using this method to avoid numpy
        # warnings from < when nans are in the array
        depth_ii = np.flatnonzero(np.isfinite(depth))  # non-nan indices
        neg_depths = np.flatnonzero(depth[depth_ii] <= 0)  # indices to depth_ii
        depth[depth_ii[neg_depths]] = np.nan

        # Remove NaN depths and truncate to when science data begins being
        # recorded and ends
        depth_ii = np.flatnonzero(np.isfinite(depth))
        # sci_indices = all_sci_indices(self.dba)  # from ooidac.processing
        # if len(sci_indices) > 0:
        #     starting_index = sci_indices[0]
        #     ending_index = sci_indices[-1]
        # else:
        #     starting_index = 0
        #     ending_index = len(timestamps) - 1
        #
        # depth_ii = depth_ii[
        #     np.logical_and(
        #         depth_ii > starting_index,
        #         depth_ii < ending_index)
        # ]

        depth = depth[depth_ii]
        ts = timestamps[depth_ii]

        if len(depth) < 2:
            logger.debug('Skipping segment that contains < 2 rows of depth')
            return

        # Create the fixed timestamp array from the min timestamp to the max
        # timestamp spaced by tsint intervals
        min_ts = ts.min()
        max_ts = ts.max()
        if max_ts - min_ts < tsint:
            logger.warning('Not enough timestamps for depth interpolation '
                           'for profile discovery')
            return

        interp_ts = np.arange(min_ts, max_ts, tsint)
        # Stretch estimated values for interpolation to span entire dataset
        interp_z = np.interp(
            interp_ts,
            ts,
            depth,
            left=depth[0],
            right=depth[-1]
        )

        filtered_z = boxcar_smooth_dataset(interp_z, filter_winsize)

        # truncate the filtered depths by science start and end times
        sci_indices = all_sci_indices(self.dba)  # from ooidac.processing
        if len(sci_indices) > 0:
            sci_start_time = timestamps[sci_indices[0]]
            sci_end_time = timestamps[sci_indices[-1]]
        else:
            sci_start_time = timestamps[0]
            sci_end_time = timestamps[-1]

        start_minimizer = abs(interp_ts - sci_start_time)
        end_minimizer = abs(interp_ts - sci_end_time)

        sci_start_ii = np.flatnonzero(
            start_minimizer == np.min(start_minimizer))
        if len(sci_start_ii) == 2:
            sci_start_ii = sci_start_ii[1]
        sci_end_ii = np.flatnonzero(end_minimizer == np.min(end_minimizer))
        if len(sci_end_ii) == 2:
            sci_end_ii = sci_end_ii[1]
        sci_truncate_indices = slice(int(sci_start_ii), int(sci_end_ii + 1))
        interp_ts = interp_ts[sci_truncate_indices]
        filtered_z = filtered_z[sci_truncate_indices]

        delta_depth = calculate_delta_depth(filtered_z)

        inflections = np.where(np.diff(delta_depth) != 0)[0] + 1
        if not inflections.any():
            return

        # inflection_times = ddts[inflections]
        inflection_times = interp_ts[inflections]  # for some reason this
        # works better
        self.inflection_times = inflection_times

        # get the first profile indices manually so that it gets all of the
        # data up to the first inflection including the inflection
        profile_i = np.flatnonzero(timestamps <= inflection_times[0])
        self._indices.append(profile_i)

        # then iterate over the inflection times and get each profile indices
        # excluding the preceding inflection and including the ending inflection
        for ii in range(1, len(inflection_times)):
            # Find all rows in the original yo that fall between the
            # interpolated timestamps
            profile_i = np.flatnonzero(
                np.logical_and(
                    timestamps > inflection_times[ii - 1],
                    timestamps <= inflection_times[ii]
                )
            )
            if len(profile_i) == 0:
                continue
            self._indices.append(profile_i)

        # lastly get the last profile manually again from the last inflection
        # time to the end of the dataset
        profile_i = np.flatnonzero(timestamps > inflection_times[-1])
        self._indices.append(profile_i)

    def find_profiles_by_depth_state(self):
        """Returns the start and stop timestamps for every profile indexed from
        the depth timeseries

        Parameters:
            time, depth

        Returns:
            A Nx2 array of the start and stop timestamps indexed from the yo

        Use filter_yo_extrema to remove invalid/incomplete profiles
        """
        profile_indexes = []

        timestamps = self.dba['llat_time']['data']

        if 'm_depth_state' not in self.dba.sensor_names:
            logging.debug('Thought there was depth state, but not')
            return profile_indexes

        depth_state = self.dba.getdata('m_depth_state')

        # remove any negative numbers or values greater than 3 since they are
        # not dive/climb/hover states and then fill the NaNs with the
        # preceeding finite value
        depth_state[np.logical_and(
            depth_state < 0.0,
            depth_state > 3.0)
        ] = np.nan
        depth_state = fwd_fill(depth_state)

        # the diff of depth state will be zero everywhere depth state stays
        # the same, but will be nonzero when states are changed.  These will
        # be the inflections.
        d_ds = np.diff(depth_state)

        inflections = np.flatnonzero(d_ds)
        inflection_times = timestamps[inflections]
        self.inflection_times = inflection_times
        # since len(timestamps) = len(d_ds)+1, the timestamps that the
        # inflection index returns will be the last timestamp for each
        # profile, so the upper times for each profile below are included by <=

        # first profile, everything up to first depth_state change. Typically
        # this is waiting on the surface to start the dive, i.e. not a
        # profile, but we eliminate that with the filters.
        profile_i = np.flatnonzero(timestamps <= inflection_times[0])
        profile_indexes.append(profile_i)

        # middle profiles
        for ii in range(1, len(inflection_times)):
            logical = np.logical_and(
                timestamps > inflection_times[ii - 1],
                timestamps <= inflection_times[ii])
            profile_i = np.flatnonzero(logical)
            profile_indexes.append(profile_i)

        # last profile
        profile_i = np.flatnonzero(timestamps > inflection_times[-1])
        profile_indexes.append(profile_i)

        # return profiled_dataset
        self._indices = profile_indexes

    def filter_profiles(self):
        """
        """
        filters = [getattr(profile_filters, func) for func in dir(
            profile_filters)
                 if func.startswith('filter')]
        profiles_to_remove = []
        for ii in range(len(self.indices)):
            profile_ii = self.indices[ii]
            profile = self.dba.slicedata(indices=profile_ii)
            for func in filters:
                remove_profile = func(profile)
                if remove_profile:
                    profiles_to_remove.append(ii)
                    logger.debug('Profile {:d} removed by {:s}'.format(
                        ii, func.__name__
                    ))
                    break
        # need to reverse profiles_to_remove indexes so that after a profile
        # is removed (popped), the remaining indices are not different for
        # the profiles_to_remove list
        profiles_to_remove.sort(reverse=True)
        for index in profiles_to_remove:
            self.indices.pop(index)


def binarize_diff(data):
    data[data <= 0] = -1
    data[data >= 0] = 1
    return data


def calculate_delta_depth(interp_data):
    delta_depth = np.diff(interp_data)
    delta_depth = binarize_diff(delta_depth)

    return delta_depth


class ProfileIndexError(IndexError):
    pass
