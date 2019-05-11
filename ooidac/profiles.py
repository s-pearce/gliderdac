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

    # ToDo: I still want to preferentially use m_depth when available for
    #  finding the profiles since there might be no science data on the
    #  climb, but I think that m_depth and the timestamps should be trimmed
    #  to start when any science data starts too, not trimmed to where
    #  science data exists, because then it would remove m_depth from climbs
    #  where there is no data.  Look at it again, but this could maybe be
    #  trimmed to the same first point as when c_profile first changes and
    #  the end trimmed to when c_battpos moves all the way forward (would
    #  need to look in reverse I think)
    def find_profiles_by_depth(self, tsint=2):
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

        timestamps = (self.dba['llat_time'] or self.dba['m_present_time'])
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
        sci_indices = all_sci_indices(self.dba)  # from ooidac.processing
        if len(sci_indices) > 0:
            starting_index = sci_indices[0]
            ending_index = sci_indices[-1]
        else:
            starting_index = 0
            ending_index = len(timestamps) - 1

        depth_ii = depth_ii[
            np.logical_and(
                depth_ii > starting_index,
                depth_ii < ending_index)
        ]

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

        filtered_z = boxcar_smooth_dataset(interp_z, int(tsint / 2))

        delta_depth = calculate_delta_depth(filtered_z)
        dts = interp_ts[:-1] + tsint / 2

        inflections = np.where(np.diff(delta_depth) != 0)[0] + 1
        if not inflections.any():
            return

        # timestamp associated with diff of delta_depth (i.e. the timestamp of
        # the 2nd derivative of depth timestamp)
        ddts = dts[:-1] + tsint / 2
        inflection_times = ddts[inflections]
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
        funcs = [getattr(profile_filters, func) for func in dir(profile_filters)
                 if func.startswith('filter')]
        profiles_to_remove = []
        for ii in range(len(self.indices)):
            profile_ii = self.indices[ii]
            profile = self.dba.slicedata(indices=profile_ii)
            for func in funcs:
                remove_profile = func(profile)
                if remove_profile:
                    profiles_to_remove.append(ii)
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
