import numpy as np
from ooidac.data_classes import DbaData
from ooidac.processing import logger
from ooidac.readers.slocum import parse_dba_header


def get_u_and_v(dba, check_files=None):
    if (
            dba.file_metadata['filename_extension'] == 'dbd'
            and check_files
            and 'm_final_water_vx' in dba.sensor_names
    ):
        vx, vy = _get_final_uv(dba, check_files)
    else:
        vx, vy = _get_initial_uv(dba)

    return vx, vy


def _get_initial_uv(dba):
    """

    :param dba:
    :return:
    """
    if 'm_initial_water_vx' in dba.sensor_names:
        vx_sensor = 'm_initial_water_vx'
        vy_sensor = 'm_initial_water_vy'
    elif 'm_water_vx' in dba.sensor_names:
        vx_sensor = 'm_water_vx'
        vy_sensor = 'm_water_vy'
    else:
        return None, None
    logger.debug('Attempting to get u & v from {:s}/vy'.format(vx_sensor))
    vx = dba[vx_sensor]
    vy = dba[vy_sensor]
    seg_initial_ii = np.flatnonzero(np.isfinite(vx['data']))
    if len(seg_initial_ii) == 0:
        vx['data'] = np.nan
        vy['data'] = np.nan
    else:
        vx['data'] = vx['data'][seg_initial_ii[-1]]
        vy['data'] = vy['data'][seg_initial_ii[-1]]

    vx.pop('sensor_name')
    vx['nc_var_name'] = 'u'
    vx['attrs']['source_sensor'] = vx_sensor
    vx['attrs']['source_file'] = dba.source_file
    vy.pop('sensor_name')
    vy['nc_var_name'] = 'v'
    vy['attrs']['source_sensor'] = vy_sensor
    vy['attrs']['source_file'] = dba.source_file

    return vx, vy


def _get_final_uv(dba, check_files):
    """return Eastward velocity `u` and Northward velocity `v` from looking
    ahead of the main glider data file into the next 2 data files given in
    the `check_files` list to retrieve `u` and `v` from the
    m_final_water_vx/vy parameter calculated 1 or 2 files/segments later.

    :param dba:
    :param check_files: sorted list of the next 2 sorted data files following
        the file being processed from the script input list.
    :return: u, v; Eastward velocity and Northward velocity in m/s as data
        particle dictionaries with metadata attributes
    """
    # ToDo: fix so that this won't fail if check_files is not

    # m_final_water_vx/vy from the current segment file applies to the
    # previous diving segment.  It is still retrieved from the current file
    # to compare when the parameter is updated in a later file.
    seg_final_vx = dba.getdata('m_final_water_vx')
    seg_final_vy = dba.getdata('m_final_water_vy')
    seg_final_ii = np.flatnonzero(np.isfinite(seg_final_vx))
    if len(seg_final_ii) == 0:
        seg_final_vx = np.nan
        seg_final_vy = np.nan
    else:
        seg_final_vx = seg_final_vx[seg_final_ii][-1]
        seg_final_vy = seg_final_vy[seg_final_ii][-1]
    mis_num = int(dba.file_metadata['the8x3_filename'][:4])
    seg_num = int(dba.file_metadata['the8x3_filename'][4:])

    # check that check_files are in the same mission and next 2 segments
    for next_dba_file in check_files:
        logger.debug(
            "Attempting to find final vx & vy in the next data file"
            "\n\t{:s}".format(next_dba_file)
        )
        header = parse_dba_header(next_dba_file)
        nxt_mis_num = int(header['the8x3_filename'][:4])
        nxt_seg_num = int(header['the8x3_filename'][4:])
        if (nxt_mis_num != mis_num or not (
                nxt_seg_num == seg_num + 1
                or nxt_seg_num == seg_num + 2)):
            logger.debug(
                'next data file {:s} is not the same mission or '
                'next 2 segments'.format(next_dba_file)
            )
            continue
        next_dba = DbaData(next_dba_file)
        if next_dba is None or next_dba.N == 0:
            continue
        if 'm_final_water_vx' not in next_dba.sensor_names:
            continue

        # get m_final_water_vx/vy from the next file
        vx = next_dba['m_final_water_vx']
        vy = next_dba['m_final_water_vy']
        vx_data = vx['data']
        vy_data = vy['data']
        next_ii = np.isfinite(vx_data)
        vx_data = vx_data[next_ii]
        vy_data = vy_data[next_ii]

        # find where the final_vx/vy value has changed from the original file
        index = np.flatnonzero(vx_data != seg_final_vx)
        if len(index) > 0:
            vx_data = vx_data[index]
            vy_data = vy_data[index]
            # if there are more than one changed values take the first one
            if len(vx_data) > 1:
                vx_data = vx_data[0]
                vy_data = vy_data[0]
            if (vx_data, vy_data) != (seg_final_vx, seg_final_vy):
                vx['data'] = vx_data
                vx.pop('sensor_name')
                vx['nc_var_name'] = 'u'
                vx['attrs']['source_sensor'] = 'm_final_water_vx'
                vx['attrs']['source_file'] = next_dba.source_file
                vy['data'] = vy_data
                vy.pop('sensor_name')
                vy['nc_var_name'] = 'v'
                vy['attrs']['source_sensor'] = 'm_final_water_vy'
                vy['attrs']['source_file'] = next_dba.source_file
                return vx, vy

    # if vx/vy not found here, return from _get_initial_uv
    return _get_initial_uv(dba)


def get_segment_time_and_pos(dba):

    # get mean segment time and nearest in time lat and lon position as
    # coordinates for the depth averaged velocities

    a = dba.ts is None
    b = dba.underwater_indices is None or len(dba.underwater_indices) == 0
    c = 'llat_latitude' not in dba.sensor_names
    d = 'llat_longitude' not in dba.sensor_names
    if a or b or c or d:
        return None, None, None

    # ToDo: dba.ts has to be m_present_time for this, should I enforce this
    #  earlier, or pair dba.timesensor with it and if it is not dba.ts,
    #  then grab m_present_time, or just grab it here regardless?
    uw_start_time = dba.ts[dba.underwater_indices[0]]
    uw_end_time = dba.ts[dba.underwater_indices[-1]]

    # borrow the sensor `llat_time`s attributes but change the data
    # to be the calculated scalar mean segment time value
    mean_segment_time = np.mean([uw_start_time, uw_end_time])
    segment_time = dba['llat_time']
    segment_time.pop('sensor_name')
    segment_time['data'] = mean_segment_time
    segment_time['nc_var_name'] = "time_uv"

    # borrow the attributes from the `llat_` sensors and replace the 'data with
    # scalar segment lat and lon
    segment_lat = dba['llat_latitude']
    segment_lat.pop('sensor_name')
    segment_lon = dba['llat_longitude']
    segment_lon.pop('sensor_name')
    lat = segment_lat['data']
    lon = segment_lon['data']

    time_diff = np.abs(dba.ts - mean_segment_time)
    closest_time_index = np.flatnonzero(time_diff == np.min(time_diff))
    seg_lat = lat[closest_time_index]
    if len(closest_time_index) > 1:
        seg_lat = np.mean(lat[closest_time_index])
        seg_lon = np.mean(lon[closest_time_index])
    else:
        seg_lat = lat[closest_time_index]
        seg_lon = lon[closest_time_index]
    segment_lat['data'] = seg_lat
    segment_lon['data'] = seg_lon
    segment_lat['nc_var_name'] = "lat_uv"
    segment_lon['nc_var_name'] = "lon_uv"

    return segment_time, segment_lat, segment_lon