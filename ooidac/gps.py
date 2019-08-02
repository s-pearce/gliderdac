import numpy as np
import logging
import os

logger = logging.getLogger(os.path.basename(__file__))


def interpolate_gps(timestamps, latitude, longitude):
    """Calculates interpolated GPS coordinates between the two surfacings
    in a single glider binary data file.

    Parameters:
        'dataset': An N by 3 numpy array of time, lat, lon pairs

    Returns interpolated gps dataset over entire time domain of dataset
    """

    dataset = np.column_stack((
        timestamps,
        latitude,
        longitude
    ))

    est_lat = np.empty((len(dataset))) * np.nan
    est_lon = np.empty((len(dataset))) * np.nan
    
    dataset = dataset[~np.isnan(dataset).any(axis=1), :]
    
    if len(dataset) == 0:
        logger.debug(
            'GPS time-series contains no valid GPS fixes for interpolation')
        return est_lat, est_lon

    # If only one GPS point, make it the same for the entire dataset
    if len(dataset) == 1:
        logger.info('Only one GPS fix, setting all records to the single fix')
        est_lat[:] = dataset[0, 1]
        est_lon[:] = dataset[0, 2]
    else:
        # Interpolate data
        est_lat = np.interp(
            timestamps,
            dataset[:, 0],
            dataset[:, 1],
            left=dataset[0, 1],
            right=dataset[-1, 1]
        )
        est_lon = np.interp(
            timestamps,
            dataset[:, 0],
            dataset[:, 2],
            left=dataset[0, 2],
            right=dataset[-1, 2]
        )

    return est_lat, est_lon


def iso2deg(iso_pos_element):
    """iso2deg converts a glider iso position element to
    a decimal degree position element.
    E.G.
    > lat = 4434.456  # a latitude, 44 deg 34.456 min
    > iso2deg(lat)
    44.574
    """
    minutes, degrees = np.modf(iso_pos_element / 100.)
    degrees = degrees + (minutes*100./60.)
    return degrees


def correct_uw_dr_pos(dba):
    """ Algorithm incomplete at this time
    A good file to test this for that has fixes in the middle of the dive is
    ce_311_2019_123_1_0.full.mrg from 311_R00007
    :param dba:
    :return:
    """
    # correct lat and lon
    gps_lon = iso2deg(dba.getdata('m_gps_lon'))
    gps_lat = iso2deg(dba.getdata('m_gps_lat'))

    dr_lat = iso2deg(dba.getdata('m_lat'))
    dr_lon = iso2deg(dba.getdata('m_lon'))
    gps_indices = np.flatnonzero(np.isfinite(gps_lon))
    gps_times = dba.ts[gps_indices]

    # this is to fix any mid segement accidental gps fixes if the glider comes
    # too close to the surface.  As these are very unlikely, i think I will have
    # to make a test data set to try this.  but I think this would work
    # ToDo: This doesn't fully work because of potential gaps during the
    #  surface times.  So I should get a segment_dive start time and a
    #  segment_dive end time and only look at gps fixes that might occur
    #  between those times for mid segment accidental fixes
    gps_dt = np.diff(gps_times)
    median_gpsdt = np.median(gps_dt)
    gaps = np.flatnonzero(gps_dt > median_gpsdt * 3)
    if len(gaps) > 1:
        logger.info("Possible mid-dive GPS fixes, attempting to remove")
        for ii in range(len(gaps) - 1):
            mid_seg_fix = gps_indices[gaps[ii] + 1:gaps[ii + 1] + 1]
            for jj in mid_seg_fix:
                if dr_lat[jj] == gps_lat[jj]:
                    lat_offset = dr_lat[jj] - dr_lat[jj - 1]  # Doesn't work
                    lon_offset = dr_lon[jj] - dr_lon[jj - 1]  # ToDo: fix this
                    next_gap_section = np.arange(gps_indices[gaps[ii + 1]] + 1,
                                                 gps_indices[gaps[ii + 1] + 1])
                    dr_lat[next_gap_section] -= lat_offset
                    dr_lon[next_gap_section] -= lon_offset
            np.delete(gps_indices, slice(gaps[ii] + 1, gaps[ii + 1] + 1))

    # now remove the remaining dr_lat/lon points that have been set to the
    # gps_fixes and any dr points that are on the surface
    dr_lat[gps_indices] = np.nan
    dr_lon[gps_indices] = np.nan
    dr_lon[dba.surface_indices] = np.nan
    dr_lat[dba.surface_indices] = np.nan

    # remove any gps indices that are underwater ( > 2 meters)
    gps_uw_ii = np.intersect1d(
        gps_indices, dba.underwater_indices,
        assume_unique=True, return_indices=True)[1]
    gps_indices = np.delete(gps_indices, gps_uw_ii)

    # remove any remaining fixes that are half a meter deeper than the median
    # fix depth
    gps_depths = dba.depth[gps_indices]
    median_gpsdepth = np.median(gps_depths)  # nans fille with interp use median
    uw_gps = np.flatnonzero(gps_depths > median_gpsdepth + 0.5)
    if len(uw_gps) > 0:
        gps_indices = np.delete(gps_indices, uw_gps)

    # all below was just pasted from the terminal, but is the basis for
    # correcting lats and lons after the above adjustments.
    pre_dive_fixes = np.intersect1d(gps_indices, dba.pre_dive_indices)
    post_dive_fixes = np.intersect1d(gps_indices, dba.post_dive_indices)

    x = [dba.ts[pre_dive_fixes[-1]], dba.ts[post_dive_fixes[0]]]
    y_lat = [gps_lat[pre_dive_fixes[-1]], gps_lat[post_dive_fixes[0]]]
    y_lon = [gps_lon[pre_dive_fixes[-1]], gps_lon[post_dive_fixes[0]]]
    p_lon = np.poly1d(np.polyfit(x, y_lon, 1))
    p_lat = np.poly1d(np.polyfit(x, y_lat, 1))

    dr_indices = np.flatnonzero(np.isfinite(dr_lon))
    dr_times = dba.ts[dr_indices]

    # method B
    new_lat = np.empty(len(dr_times))
    new_lon = np.empty(len(dr_times))
    initial_new_lat = p_lat(dr_times[0])
    initial_new_lon = p_lon(dr_times[0])
    final_new_lat = p_lat(dr_times[-1])
    final_new_lon = p_lon(dr_times[-1])
    new_lat[0] = initial_new_lat
    new_lon[0] = initial_new_lon
    new_lat[-1] = final_new_lat
    new_lon[-1] = final_new_lon

    init_offset_lat = dr_lat[dr_indices[0]] - initial_new_lat
    init_offset_lon = dr_lon[dr_indices[0]] - initial_new_lon
    final_offset_lat = dr_lat[dr_indices[-1]] - final_new_lat
    final_offset_lon = dr_lon[dr_indices[-1]] - final_new_lon

    p_offset_lat = np.poly1d(np.polyfit(
        [dr_times[0], dr_times[-1]],
        [init_offset_lat, final_offset_lat],
        1))
    p_offset_lon = np.poly1d(np.polyfit(
        [dr_times[0], dr_times[-1]],
        [init_offset_lon, final_offset_lon],
        1))

    lat_offsets = p_offset_lat(dr_times[1:-1])
    lon_offsets = p_offset_lon(dr_times[1:-1])

    new_lat[1:-1] = dr_lat[dr_indices[1:-1]] - lat_offsets
    new_lon[1:-1] = dr_lon[dr_indices[1:-1]] - lon_offsets

    lats = np.empty(len(dba.ts))
    lons = np.empty(len(dba.ts))
    lats[dr_indices] = new_lat
    lons[dr_indices] = new_lon
    lats[gps_indices] = gps_lat[gps_indices]
    lons[gps_indices] = gps_lon[gps_indices]

    positions_ii = np.union1d(gps_indices, dr_indices)
    lats = np.interp(
        dba.ts, dba.ts[positions_ii], lats[positions_ii],
        left=lats[positions_ii[0]], right=lats[positions_ii[-1]])
    lons = np.interp(
        dba.ts, dba.ts[positions_ii], lons[positions_ii],
        left=lons[positions_ii[0]], right=lons[positions_ii[-1]])
    return lats, lons
