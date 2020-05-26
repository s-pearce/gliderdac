import os
import logging
import numpy as np

from configuration import DATA_CONFIG_LIST, REQUIRED_SENSORS
from configuration import MIN_DATA_VALS, MIN_DIVE_DEPTH
from configuration import DAV_SENSORS, DEPTH_SENSOR
from ooidac.constants import SLOCUM_SALINITY_SENSORS

logger = logging.getLogger(os.path.basename(__name__))


class FileCheck(object):
    required_sensors = False
    any_science_data = False
    avail_sci_data = []
    dav_sensors = False
    file_good = False


def check_file_goodness(gldata):
    """The main goodness check function.
    Checks for required sensors, available science sensors, """
    fc = FileCheck()
    fc.required_sensors = check_required_sensors(gldata)
    fc.avail_sci_data = sci_data_available(gldata)
    fc.dav_sensors = check_for_dav_sensors(gldata)
    if len(fc.avail_sci_data) > 0:
        fc.any_science_data = True
    fc.file_good = fc.required_sensors and fc.any_science_data

    return fc


def check_required_sensors(gldata):
    required_sensors_present = True
    for sensor in REQUIRED_SENSORS:
        if sensor not in gldata.sensor_names:
            required_sensors_present = False
            logger.warning('Required Sensor: {:s} not present in {:s}'.format(
                sensor, gldata.source_file)
            )
    return required_sensors_present


def check_for_dav_sensors(gldata):
    sensor_names = gldata.sensor_names
    # dav_sensors = [  # replaced by explicit strings in configuration.py
    #     ('m_final_water_vx', 'm_final_water_vy'),
    #     ('m_water_vx', 'm_water_vy'),
    #     ('m_initial_water_vx', 'm_initial_water_vy')
    # ]
    check = []
    for vx, vy in DAV_SENSORS:
        if vx in sensor_names and vy in sensor_names:
            check.append((vx, vy))
    if len(check) > 0:
        dav_exists = True
    else:
        dav_exists = False
    return dav_exists, check


def check_if_dive(gldata):
    diving_segment = False
    depth = gldata.getdata(DEPTH_SENSOR)
    max_depth = np.nanmax(depth)
    if max_depth > MIN_DIVE_DEPTH:
        diving_segment = True
    return diving_segment


def check_for_any_sci_data(gldata):
    data_exists = False
    for sensor in DATA_CONFIG_LIST:
        if sensor in gldata.sensor_names:
            data = gldata.getdata(sensor)
            if len(np.flatnonzero(np.isfinite(data))) > MIN_DATA_VALS:
                data_exists = True
                break
    return data_exists


def sci_data_available(gldata):
    sci_sensors = []
    for sensor in DATA_CONFIG_LIST:
        if sensor in gldata.sensor_names:
            data = gldata.getdata(sensor)
            if len(np.flatnonzero(np.isfinite(data))) > MIN_DATA_VALS:
                # config
                sci_sensors.append(sensor)
            else:
                logger.debug(
                    'Science data {:s} has less than {:d} values in file '
                    '{:s}'.format(sensor, MIN_DATA_VALS, gldata.source_file)
                )
        else:
            # ToDo: Would be good to eventually get the acutal configuration
            #       file path here for the eventual multiple configuration files
            logger.warning(
                'Science Sensor {:s} not found in data file {:s}, but is '
                'found in configuration.py'.format(sensor, gldata.source_file)
            )
    return sci_sensors


# The function below is deprecated by `check_required_sensors` above
# def check_ctd_sensors(gldata):
#     ctd_sensors_exist = False
#     ctd_sensors = np.intersect1d(
#         DATA_CONFIG_LIST, [
#             'sci_water_cond', 'sci_water_temp',
#             'm_water_cond', 'm_water_temp']
#     )
#     for sensor in ctd_sensors:
#         if sensor not in gldata.sensor_names:
#             logging.warning(
#                 ('Sensor {:s} for processing CTD data not found in '
#                  'data file {:s}').format(sensor, gldata.source_file)
#             )
#             return False
#     return True


def cond_sensor(gldata):
    for sensor in SLOCUM_SALINITY_SENSORS:
        if sensor in gldata.sensor_names:
            return sensor
    return None
