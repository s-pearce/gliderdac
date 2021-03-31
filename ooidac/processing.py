
import os
import logging
from copy import deepcopy
import gsw
import ooidac.gps as gps
import numpy as np
from ooidac.data_classes import DbaData
from ooidac.readers.slocum import parse_dba_header
from ooidac.utilities import fwd_fill
from ooidac.ctd import calculate_practical_salinity, calculate_density
from ooidac.processing_dir.fluorometer import flo_bback_total, flo_beta
from ooidac.processing_dir.oxygen_calculation import calc_o2, do2_SVU
from configuration import SCITIMESENSOR
from ooidac.constants import (
    SLOCUM_TIMESTAMP_SENSORS,
    SLOCUM_PRESSURE_SENSORS,
    SLOCUM_DEPTH_SENSORS)
from configuration import DATA_CONFIG_LIST

logger = logging.getLogger(os.path.basename(__name__))


def create_llat_sensors(
        dba, timesensor=None, pressuresensor=None,
        depthsensor=None, z_from_p=True):

    # List of available dba sensors
    dba_sensors = dba.sensor_names

    # Select the time sensor
    time_sensor = select_time_sensor(dba, timesensor=timesensor)
    if not time_sensor:
        return
    # Select the pressure sensor
    pressure_sensor = select_pressure_sensor(dba, pressuresensor=pressuresensor)
    # Select the depth sensor
    depth_sensor = select_depth_sensor(dba, depthsensor=depthsensor)
    # We must have either a pressure_sensor or depth_sensor to continue
    if not pressure_sensor and not depth_sensor:
        logger.warning(
            'No pressure sensor and no depth sensor found: {:s}'.format(
                dba.source_file)
        )
        return

    # Must have m_gps_lat and m_gps_lon to convert to decimal degrees
    if 'm_gps_lat' not in dba_sensors or 'm_gps_lon' not in dba_sensors:
        logger.warning(
            'Missing m_gps_lat and/or m_gps_lon: {:s}'.format(dba.source_file)
        )
        return

    lat_sensor, lon_sensor = lat_and_lon_coordinates(dba, time_sensor)

    # If no depth_sensor was selected, use llat_latitude, llat_longitude
    # and llat_pressure to calculate
    if not depth_sensor or z_from_p:
        if pressure_sensor:
            logger.debug(
                'Calculating depth from selected pressure sensor: {:s}'.format(
                    pressure_sensor['attrs']['source_sensor']))

            depth_sensor = {
                'sensor_name': 'llat_depth',
                'attrs': {
                    'source_sensor': 'llat_pressure,llat_latitude',
                    'comment': (
                        u'Calculated from llat_pressure and '
                        u'llat_latitude using gsw.z_from_p'
                    )
                },
                'data': -gsw.z_from_p(
                    pressure_sensor['data'], lat_sensor['data'])
                }
        else:
            logging.warning(
                'No pressure sensor found for calculating depth')

    # Append the llat variables
    dba.add_data(time_sensor)

    if pressure_sensor:
        dba.add_data(pressure_sensor)

    dba.add_data(depth_sensor)

    dba.add_data(lat_sensor)

    dba.add_data(lon_sensor)

    return dba


# ToDo: parts or all of this would be better suited in a configuration
#  class/module
def select_time_sensor(dba, timesensor=None):
    # Figure out which time sensor to select
    if timesensor:
        if timesensor not in dba.sensor_names:
            logger.warning(
                'Specified timesensor {:s} not found in dba: {:s}, '
                'auto-choosing one instead'.format(
                    timesensor, dba.file_metadata['source_file']))
            timesensor = _autochoose(dba, SLOCUM_TIMESTAMP_SENSORS, 'time')
    else:
        timesensor = _autochoose(dba, SLOCUM_TIMESTAMP_SENSORS, 'time')

    if not timesensor:
        return

    time_sensor = deepcopy(dba[timesensor])
    if not time_sensor:
        return
    time_sensor['sensor_name'] = 'llat_time'
    time_sensor['attrs']['source_sensor'] = timesensor
    time_sensor['attrs']['comment'] = u'Alias for {:s}'.format(timesensor)
    time_sensor['attrs']['units'] = 'seconds since 1970-01-01 00:00:00Z'

    # the new sensor data is the same as the source sensor data, so no ['data']
    # update needed

    return time_sensor


# ToDo: parts or all of this would be better suited in a configuration class
#  or module
def select_pressure_sensor(dba, pressuresensor=None):
    """Returns selected pressure sensor name and pressure array in decibars"""

    # List of available dba sensors
    dba_sensors = dba.sensor_names

    # User pressuresensor if specified
    if pressuresensor:
        if pressuresensor not in dba_sensors:
            logger.warning(
                'Specified pressuresensor {:s} not found in dba {:s}. '
                'Auto-choosing pressure sensor instead'.format(
                    pressuresensor, dba.file_metadata['source_file']))
            pressuresensor = _autochoose(
                dba, SLOCUM_PRESSURE_SENSORS, 'pressure'
            )
    else:
        pressuresensor = _autochoose(dba, SLOCUM_PRESSURE_SENSORS, 'pressure')

    if not pressuresensor:
        return

    pressure_sensor = deepcopy(dba[pressuresensor])
    if not pressure_sensor:
        return
    pressure_sensor['sensor_name'] = 'llat_pressure'
    pressure_sensor['attrs']['source_sensor'] = pressuresensor
    pressure_sensor['attrs']['comment'] = (
        u'Alias for {:s}, multiplied by 10 to convert from bar to dbar'.format(
            pressuresensor)
    )

    # Convert the pressure sensor from bar to dbar
    pressure_sensor['data'] = pressure_sensor['data'] * 10
    pressure_sensor['attrs']['units'] = 'dbar'

    return pressure_sensor


# ToDo: Parts or all of this would be better suited to a configuration class
#  or module
def select_depth_sensor(dba, depthsensor=None):
    # List of available dba sensors
    dba_sensors = dba.sensor_names

    # User pressuresensor if specified
    if depthsensor:
        if depthsensor not in dba_sensors:
            logger.warning(
                'Specified depthsensor {:s} not found in dba: {:s}, '
                'auto-choosing depth sensor instead'.format(
                    depthsensor, dba.file_metadata['source_file'])
            )
            depthsensor = _autochoose(dba, SLOCUM_DEPTH_SENSORS, 'depth')
    else:
        depthsensor = _autochoose(dba, SLOCUM_DEPTH_SENSORS, 'depth')

    if not depthsensor:
        return

    depth_sensor = deepcopy(dba[depthsensor])
    if not depth_sensor:
        return
    depth_sensor['sensor_name'] = 'llat_depth'
    depth_sensor['attrs']['source_sensor'] = depthsensor
    depth_sensor['attrs']['comment'] = u'Alias for {:s}'.format(depthsensor)

    # the new sensor data is the same as the source sensor data, so no ['data']
    # update needed

    return depth_sensor


def lat_and_lon_coordinates(dba, time_sensor):
    # Convert m_gps_lat to decimal degrees and create the new sensor
    # definition
    lat_sensor = deepcopy(dba['m_gps_lat'])
    lat_sensor['sensor_name'] = 'llat_latitude'
    lat_sensor['attrs']['source_sensor'] = u'm_gps_lat'

    # Skip default values (69696969)
    # ToDo: fix this so it doesn't print a warning
    lats = lat_sensor['data']
    lat_ii = np.flatnonzero(np.isfinite(lats))
    bad_lats = lats[lat_ii] > 9000.0
    lats[lat_ii[bad_lats]] = np.nan
    lat_sensor['data'] = gps.iso2deg(lats)

    # lat_sensor['data'][lat_sensor['data'] > 9000.0] = np.nan
    # lat_sensor['data'] = gps.iso2deg(lat_sensor['data'])

    # Convert m_gps_lon to decimal degrees and create the new sensor
    # definition
    lon_sensor = deepcopy(dba['m_gps_lon'])
    lon_sensor['sensor_name'] = 'llat_longitude'
    lon_sensor['attrs']['source_sensor'] = u'm_gps_lon'

    # Skip default values (69696969)
    # ToDo: fix this so it doesn't print a warning
    lons = lon_sensor['data']
    lon_ii = np.flatnonzero(np.isfinite(lons))
    bad_lons = lons[lon_ii] > 18000.0
    lons[lon_ii[bad_lons]] = np.nan
    lon_sensor['data'] = gps.iso2deg(lons)

    # lon_sensor['data'][lon_sensor['data'] > 18000] = np.nan
    # lon_sensor['data'] = gps.iso2deg(lon_sensor['data'])

    logging.info('Filling lat and lon coordinates by interpolation '
                 'between GPS fixes')
    # Interpolate llat_latitude and llat_longitude
    lat_sensor['data'], lon_sensor['data'] = gps.interpolate_gps(
        time_sensor['data'], lat_sensor['data'], lon_sensor['data']
    )
    lat_sensor['attrs']['comment'] = (
        u'm_gps_lat converted to decimal degrees and interpolated'
    )
    lon_sensor['attrs']['comment'] = (
        u'm_gps_lon converted to decimal degrees and interpolated'
    )

    return lat_sensor, lon_sensor


def _autochoose(dba, sensorlist, sensortype):
    sensor = None
    for s in sensorlist:
        if s in dba.sensor_names:
            sensor = s
            logger.info('Auto-chose {:s} sensor: {:s}'.format(
                sensortype, sensor)
            )
            break
    if not sensor:
        logger.warning(
            'No {:s} sensor found in dba: {:s}'.format(
                sensortype, dba.file_metadata['source_file'])
        )
    return sensor


def pitch_and_roll(dba, fill='fwd fill'):
    """adds new sensors `pitch` and `roll` to a GliderData instance from
    `m_pitch` and `m_roll` converted to degrees from radians.

    :param dba:  A GliderData or DbaData instance
    :param fill:
    :return: dba:  The same GliderData instance with `pitch` and `roll` added
    """
    if fill == 'fwd fill':
        fill_function = fwd_fill
    elif fill == 'interp':
        def fill_function(param):
            filled_param = np.interp(
                dba.ts, dba.ts[np.isfinite(param)], param[np.isfinite(param)])
            return filled_param
    else:
        def fill_function(param):
            return param

    pitch = dba['m_pitch']
    pitch['sensor_name'] = 'pitch'
    pitch['attrs']['units'] = 'degrees'
    pitch['attrs']['comment'] = (
        'm_pitch converted to degrees and forward filled')
    pitch['data'] = fill_function(np.degrees(pitch['data']))

    roll = dba['m_roll']
    roll['sensor_name'] = 'roll'
    roll['attrs']['units'] = 'degrees'
    roll['attrs']['comment'] = (
        'm_roll converted to degrees and forward filled')
    roll['data'] = fwd_fill(np.degrees(roll['data']))

    dba.add_data(pitch)
    dba.add_data(roll)

    return dba


# TODO: Build a sensor defs (or separate attributes classes) that clearly
#  read in the JSON definitions files and can use those rather than a murky
#  singular ncw class.  Then that can be passed to this function.
def ctd_data(dba, ctd_sensors):
    # before beginning, check that all the proper sensors are there
    for sensor in ctd_sensors:
        if sensor not in dba.sensor_names:
            logging.warning(
                ('Sensor {:s} for processing CTD data not found in '
                 'dba file {:s}').format(sensor, dba.source_file)
            )
            return

    # if that didn't return, get the sensors needed
    pres = dba['llat_pressure']
    lat = dba['llat_latitude']
    lon = dba['llat_longitude']
    temp = dba['sci_water_temp'] or dba['m_water_temp']
    cond = dba['sci_water_cond'] or dba['m_water_cond']

    # make sure none of the variables are completely empty of data
    for var in [pres, lat, lon, temp, cond]:
        if np.all(np.isnan(var['data'])):
            logging.warning(
                'dba file {:s} contains no valid {:s} values'.format(
                    dba.source_file,
                    var['sensor_name'])
            )
            return

    # Calculate mean llat_latitude and mean llat_longitude
    mean_lat = np.nanmean(lat['data'])
    mean_lon = np.nanmean(lon['data'])

    # Calculate practical salinity
    prac_sal = calculate_practical_salinity(
        cond['data'], temp['data'], pres['data'])
    # Add salinity to the dba
    dba['salinity'] = {
        'sensor_name': 'salinity',
        'attrs': {},  # these get filled in later by the netCDF writer
        'data': prac_sal
    }

    # Calculate density
    density = calculate_density(
        temp['data'], pres['data'], prac_sal, mean_lat, mean_lon)
    # Add density to the dba
    dba['density'] = {
        'sensor_name': 'density',
        'attrs': {},  # these get filled in later by the netCDF writer
        'data': density
    }
    return dba


class CTDprocessingError(Exception):
    pass


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


# ToDo: make these return an array instead?


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


def backscatter_total(gldata, src_sensor, var_name, **kwargs):
    """Calculate total backscatter

    Can either calculate total backscatter from the scaled units uE/m^2sec or
    from raw instrument counts if `dark_counts` and `scale_factor` are
    included as inputs.  If inputs `wlngth`, `theta`, and `xfactor` are not
    included, the parameters for a 3 channel ECO triplet with 700 nm
    wavelength are assumed (e.g. FLBBCD)

    :param gldata: a GliderData instance
    :param str src_sensor: The source sensor(variable) to use for the
        calculation
    :param str var_name: The variable name to use for backscatter added to
        the GliderData instance.
    :param kwargs: optional keyword arguments
        wlngth: Wavelength of the backscatter sensor [nm]. Default is 700 nm
        theta: Centroid Angle of the backscatter sensor [degrees].
            Default is 124 degrees.
        xfactor: Chi-factor of the backscatter sensor.  Default is 1.076
        dark_counts : The raw counts value for the sensor in the dark from the
            calibration information. Must be included with `scale_factor`,
            otherwise ignored
        scale_factor: The wet scale factor for the sensor from the calibration
            information for the appropriate end units. Must be included with
            `dark_counts` otherwise ignored

    :return: The GliderData instance with the backscatter variable added.
    """
    # set defaults if these arguments are not included in the kwargs
    required_args = ['theta', 'wlngth', 'xfactor']
    if (
            'wlngth' not in kwargs and
            'theta' not in kwargs and
            'xfactor' not in kwargs):
        kwargs['wlngth'] = 700.0
        kwargs['theta'] = 124.0
        kwargs['xfactor'] = 1.076
    assert all([x in kwargs for x in required_args]), (
        '"wlngth", "theta", and "xfactor" must either all be present or absent '
        'in the "processing" dictionary for "{:s}".'.format(var_name)
    )

    # return if sensor is not found in GliderData instance
    if src_sensor not in gldata.sensor_names:
        logger.warning(
            'Adding {:s} as all NaNs because source variable {:s} for '
            'processing does not exist'.format(
                var_name, src_sensor))
        gldata = _add_nan_variable(gldata, var_name)
        return gldata

    # copy the particle structure to keep the metadata from the src_sensor
    backscatter_particle = deepcopy(gldata[src_sensor])

    # if src_sensor ends with _sig it is in raw counts and needs to be converted
    # to beta using the calibration values for dark_counts and scale_factor,
    # which both need to be present.
    if src_sensor.endswith('_sig'):
        assert ('dark_counts' in kwargs and 'scale_factor' in kwargs), (
            '"dark_counts" and "scale_factor" must be included as "processing" '
            'dictionary arguments for "{:s}" in sensor_defs.json.'.format(
                var_name)
        )
        dark_counts = kwargs.pop('dark_counts')
        scale_factor = kwargs.pop('scale_factor')
        counts = backscatter_particle['data'].copy()
        beta = flo_beta(counts, dark_counts, scale_factor)
    else:
        beta = backscatter_particle['data'].copy()

    backscatter_particle['data'] = np.full(len(beta), np.nan)

    beta_ii = np.isfinite(beta)
    timestamps = gldata.getdata('m_present_time')
    beta_ts = timestamps[beta_ii]
    beta = beta[beta_ii]
    temp = gldata.getdata('sci_water_temp')
    salt = gldata.getdata('salinity')

    temp_bt = np.interp(
        beta_ts, timestamps[np.isfinite(temp)], temp[np.isfinite(temp)])
    salt_bt = np.interp(
        beta_ts, timestamps[np.isfinite(salt)], salt[np.isfinite(salt)])

    theta = kwargs['theta']
    wlngth = kwargs['wlngth']
    xfactor = kwargs['xfactor']
    bback = flo_bback_total(beta, temp_bt, salt_bt, theta, wlngth, xfactor)

    backscatter_particle['data'][beta_ii] = bback
    backscatter_particle['sensor_name'] = var_name
    backscatter_particle['attrs']['units'] = 'm-1'

    gldata.add_data(backscatter_particle)
    return gldata


def reduce_to_sci_data(gldata):
    """Reduces a GliderData instance to timestamps that only exhibit at least
    one measurement of science instrument data as determined by the
    configuration parameter DATA_CONFIG_LIST

    :param gldata: A GliderData instance
    :return: A new GliderData instance reduced to only science data
    """
    sci_indices = all_sci_indices(gldata)

    reduced_gldata = gldata.slicedata(indices=sci_indices)
    return reduced_gldata


def all_sci_indices(gldata):
    """

    :param gldata: A GliderData instance
    :return:
    """
    sci_indices = np.array([], dtype=np.int64)
    for sci_sensor in DATA_CONFIG_LIST:
        sci_data = gldata.getdata(sci_sensor)
        sci_ii = np.flatnonzero(np.isfinite(sci_data))
        sci_indices = np.union1d(sci_indices, sci_ii)
    return sci_indices


def remove_sci_init_zeros(gldata, available_sensors):
    """Sets values to NaNs where all primary science data sensors equal 0.0 at
    a timestamp, such as might occur when a science bay initializes.
    `available_sensors` is a list of available science data sensors.  Typically
    these would be returned from calling the `check_file_goodness` function from
    the `data_checks.py` module which compares available sensors to the
    DATA_CONFIG_LIST in the configuration file.

    :param gldata: GliderData instance
    :param available_sensors: list of available glider sensors based on the
        check_file_goodness test in data_checks.py
    :return: The same GliderData instance with any initialization zeros
    changed to NaNs
    """
    zeros_ii = np.array([])
    for sensor in available_sensors:
        var = gldata.getdata(sensor)
        if sensor == available_sensors[0]:
            # get the zeros from the first sensor differently to start the
            # set intersection below
            zeros_ii = np.flatnonzero(var == 0.0)
        else:
            var_zero_ii = np.flatnonzero(var == 0.0)
            # intersection will find the timestamp where all of the science
            # sensors are zero
            zeros_ii = np.intersect1d(zeros_ii, var_zero_ii)
    if len(zeros_ii) > 0:
        gldata.update_data(available_sensors, zeros_ii, np.nan)
    return gldata


def recalc_chlor(dba, dark_offset, scale_factor):
    """ Recalculates chlorophyll from the raw signal

    :param dba: GliderData instance
    :param dark_offset: Dark Offset calibration parameter
    :param scale_factor: Scale Factor calibration parameter
    :return: The GliderData instance with the new parameter added
    """
    if 'sci_flbbcd_chlor_sig' not in dba.sensor_names:
        logger.warning(
            "sci_flbbcd_chlor_sig is not present to recalculate chlorophyll. "
            "Filling with NaNs instead.")
        # This needs to return a variable called `corrected_chlor` or else the
        # code will not continue with the rest of the potentially good data, so
        # this should return `corrected_chlor` full of nans since it can not be
        # corrected
        nanchlor = deepcopy(dba['sci_flbbcd_chlor_units'])
        nanchlor['sensor_name'] = "corrected_chlor"
        nanchlor['data'] = np.full(len(dba), np.nan)
        dba.add_data(nanchlor)
        return _add_nan_variable(
                dba, "corrected_chlor", "sci_flbbcd_chlor_units")
    chlor_sig = dba.getdata('sci_flbbcd_chlor_sig')
    chlor_units = deepcopy(dba['sci_flbbcd_chlor_units'])
    new_chlor = scale_factor * (chlor_sig - dark_offset)
    chlor_units['data'] = new_chlor
    chlor_units['attrs']['comment'] = (
        "Chlorophyll recalculated from signal using calibration parameters")
    chlor_units['sensor_name'] = "corrected_chlor"
    chlor_units['attrs']['source_sensor'] = "sci_flbbcd_chlor_sig"
    dba.add_data(chlor_units)

    return dba


def recalc_par(dba, sensor_dark, scale_factor):
    """ Recalculates PAR from the raw signal

    :param dba: GliderData instance
    :param sensor_dark: Dark Offset calibration parameter from the calibration
        certificate
    :param scale_factor: Appropriate units Wet Scale Factor calibration
        parameter from the calibration certificate
    :return: The GliderData instance with the new parameter added
    """
    if 'sci_bsipar_sensor_volts' not in dba.sensor_names:
        logger.warning(
            "sci_bsipar_sensor_volts is not present to recalculate PAR. "
            "Filling with NaNs instead.")
        # This needs to return a variable called `corrected_par` or else the
        # code will not continue with the rest of the potentially good data, so
        # this should return `corrected_par` full of nans since it can not be
        # corrected
        return _add_nan_variable(dba, "corrected_par", "sci_bsipar_par")
    par_volts = dba.getdata('sci_bsipar_sensor_volts')
    # remove the initialization where sensor volts == 0.0
    par_volts[par_volts == 0.0] = np.nan
    par_units = deepcopy(dba['sci_bsipar_par'])
    new_par = (par_volts - sensor_dark) / scale_factor
    par_units['data'] = new_par
    par_units['attrs']['comment'] = (
        "PAR recalculated from signal using calibration parameters")
    par_units['sensor_name'] = "corrected_par"
    par_units['attrs']['source_sensor'] = "sci_bsipar_sensor_volts"
    dba.add_data(par_units)

    return dba


def check_and_recalc_o2(dba, calc_type, cal_dict):
    """ Checks the oxygen calculation and re-calculates from calphase if
    salinity was entered as 35 instead of 0 and uses the real gas constant if
    the ideal gas constant was used for MkII calculations in firmware earlier
    than 4.5.7.

    :param dba: GliderData instance
    :param calc_type
    :param cal_dict:
    :return: The GliderData instance with the new parameter added
    """
    if 'sci_oxy4_calphase' not in dba.sensor_names:
        logger.warning(
            "sci_oxy4_calphase is not present to recalculate oxygen. Filling "
            "with NaNs instead.")
        # This needs to return a variable called `corrected_oxygen` or else the
        # code will not continue with the rest of the potentially good data, so
        # this should return `corrected_oxygen` full of nans since it can not be
        # corrected
        return _add_nan_variable(dba, 'corrected_oxygen', 'sci_oxy4_oxygen')
    calphase = dba.getdata('sci_oxy4_calphase')
    oxytemp = dba.getdata('sci_oxy4_temp')
    bads = calphase == 0.0
    calphase[bads] = np.nan
    oxytemp[bads] = np.nan
    if calc_type == 'SVU':
        csv = cal_dict['SVUFoilCoef']
        conc_coef = cal_dict['ConcCoef']
        new_oxy = do2_SVU(calphase, oxytemp, csv, conc_coef, salt=0.0)
    elif calc_type == 'MkII':
        C = cal_dict['C']
        fpt = cal_dict['FoilPolyDegT']
        fpo = cal_dict['FoilPolyDegO']
        conc_coef = cal_dict['ConcCoef']
        new_oxy = calc_o2(calphase, oxytemp, C, M=fpt, N=fpo, cc=conc_coef)[0]
    else:
        logger.warning(
            "No oxygen calculation type configured for re-calculation.  Add "
            "'calculation_type': 'SVU' or 'MkII' to the 'corrected_oxygen' "
            "section in sensor_defs.json")
        dba = _add_nan_variable(dba, "corrected_oxygen", "sci_oxy4_oxygen")
        return dba

    oxy_units = deepcopy(dba['sci_oxy4_oxygen'])
    oxy_units['data'] = new_oxy
    oxy_units['attrs']['comment'] = (
        "Oxygen recalculated from signal using calibration parameters")
    oxy_units['sensor_name'] = "temp_corrected_oxygen"
    oxy_units['attrs']['source_sensor'] = "sci_oxy4_calphase, sci_oxy4_temp"
    dba.add_data(oxy_units)

    return dba


def o2_s_and_p_comp(dba, o2sensor='sci_oxy4_oxygen'):
    if o2sensor not in dba.sensor_names:
        logger.warning(
            'Oxygen data not found in data file {:s}'.format(dba.source_file)
        )
        return dba

    oxygen = dba[o2sensor]
    oxy = oxygen['data'].copy()
    timestamps = dba.getdata(SCITIMESENSOR)
    sp = dba.getdata('salinity')
    p = dba.getdata('llat_pressure')
    t = dba.getdata('sci_water_temp')

    oxy_ii = np.isfinite(oxy)
    oxy = oxy[oxy_ii]
    oxy_ts = timestamps[oxy_ii]

    sp = np.interp(oxy_ts, timestamps[np.isfinite(sp)], sp[np.isfinite(sp)])
    p = np.interp(oxy_ts, timestamps[np.isfinite(p)], p[np.isfinite(p)])
    t = np.interp(oxy_ts, timestamps[np.isfinite(t)], t[np.isfinite(t)])

    lon = dba.getdata('llat_longitude')[oxy_ii]  # should already be interp'ed
    lat = dba.getdata('llat_latitude')[oxy_ii]

    # density calculation from GSW toolbox
    sa = gsw.SA_from_SP(sp, p, lon, lat)
    ct = gsw.CT_from_t(sa, t, p)
    pdens = gsw.rho(sa, ct, 0.0)  # potential referenced to p=0

    # Convert from volume to mass units:
    do = 1000*oxy/pdens

    # Pressure correction:
    do = (1 + (0.032*p)/1000) * do

    # Salinity correction (Garcia and Gordon, 1992, combined fit):
    s0 = 0
    ts = np.log((298.15-t)/(273.15+t))
    b0 = -6.24097e-3
    b1 = -6.93498e-3
    b2 = -6.90358e-3
    b3 = -4.29155e-3
    c0 = -3.11680e-7
    bts = b0 + b1*ts + b2*ts**2 + b3*ts**3
    do = np.exp((sp-s0)*bts + c0*(sp**2-s0**2)) * do

    oxygen['sensor_name'] = 'oxygen'
    oxygen['data'] = np.full(len(oxy_ii), np.nan)
    oxygen['data'][oxy_ii] = do
    oxygen['attrs']['units'] = "umol kg-1"
    if 'comment' in oxygen['attrs']:
        comment = oxygen['attrs']['comment'] + "; "
    else:
        comment = ''
    oxygen['attrs']['comment'] = comment + (
        "Oxygen concentration has been compensated for salinity and "
        "pressure, but has not been corrected for the depth offset "
        "due to pitch of the glider and sensor offset from the CTD.")
    dba.add_data(oxygen)

    return dba


def replace_missing_sensors(gldata, available_sensors):
    """Replaces missing variables/sensors from the gliderdata object so things
    don't break if a variable is missing.

    This will compare the `available_sensors` input (which should be a subset of
    DATA_CONFIG_LIST from the gdac configuration file) with DATA_CONFIG_LIST,
    and if any sensors are missing, it will create it as an array of NaNs,
    just so the code does not break.  This is just considered a patch until
    better code exists.

    :param gldata: GliderData instance
    :param available_sensors: a list of the sensors available as a subset of
    the DATA_CONFIG_LIST from the gdac configuration file
    :return: The GliderData instance input with the new sensor added
    """
    for sensor in DATA_CONFIG_LIST:
        if sensor not in available_sensors:
            new_sensor = {
                "sensor_name": sensor,
                "data": np.full(len(gldata), np.nan),
                "attrs": {
                    "units": "N/A",
                    "bytes": 4,
                    "comment": (
                        "Data variable {:s} is absent from this data "
                        "file and is all _FillValue values".format(
                            sensor)),
                    "source_sensor": sensor,
                    "long_name": sensor}
            }
            gldata.add_data(new_sensor)
    return gldata


def _add_nan_variable(gdata, variable_name, derived_variable=None):
    """ Adds a variable named `variable_name` to the glider data object
    `gdata` that is filled with all NaNs using the attributes from variable /
    sensor `derived_variable`

    :param gdata: GliderData instance to add the NaN variable to.
    :param variable_name: string
        The name of the new NaN-filled variable to add.
    :param derived_variable: string
        The name of the variable to use to derive the attributes from such as
        units, long_name, and bytes (that indicates data type).
    :return: gdata: GliderData instance with new NaN-filled variable added
    """
    if derived_variable is not None and derived_variable in gdata.sensor_names:
        nanvar = deepcopy(gdata[derived_variable])
    else:
        nanvar = {'attrs': {}}
    nanvar['sensor_name'] = variable_name
    nanvar['data'] = np.full(len(gdata), np.nan)
    gdata.add_data(nanvar)
    return gdata
