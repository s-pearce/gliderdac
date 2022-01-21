
import os
import logging
from copy import deepcopy
import gsw
import numpy as np
from ooidac.processing.gps import lat_and_lon_coordinates
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


def init_processing_func(func_dict):
    """Import a generic processing function

    """
    module = func_dict['module']

    # check if the module has path elements in it and separate those out.
    # Check the paths exist if present.
    if "\\" in module or "/" in module:
        module_path = os.path.abspath(module)
        assert os.path.exists(module_path), (
            "Processing module {:s} not found".format(module_path))
        module_dir = os.path.dirname(module_path)
    elif 'module_dir' in func_dict:
        module_dir = func_dict['module_dir']
        assert os.path.isdir(module_dir), (
            "Processing module directory {:s} not found".format(module_dir))
    else:
        module_dir = None

    # if there were path elements in `module` and add the path to sys.path
    # for importing if not already there
    if module_dir and module_dir not in sys.path:
        sys.path.append(module_dir)

    # if there were path elements, these will remove them and do nothing if not.
    module = os.path.basename(module)
    modulename = module.replace('.py', '')

    # import the module
    try:
        mod = importlib.import_module(modulename)
    except ModuleNotFoundError:
        logger.error(
            "Processing Module from user input {:s} not found. Please "
            "re-edit the processing step and try running "
            "again.".format(modulename))
        sys.exit(1)

    # get the specific processing function to keep in the processing dictionary
    func_name = func_dict['function']
    try:
        func = getattr(mod, func_name)
    except AttributeError:
        logger.error(
            "Function {:s} not found in module {:s}. Please edit the processing"
            " info and try running again.".format(func_name, modulename)
        )
        sys.exit()

    func_dict['func'] = func
    return func_dict


def init_processing_dict(var_processing_dict):
    """Initialize the processing dictionary by importing the generic functions
    """

    # pre-emptively insert the processing directory into the module search path
    sys.path.insert(0, PROCESSING_DIR)

    # each processing instruction block, import the function and store the
    # function as an object in the processing dictionary
    for var, func_dict in var_processing_dict.items():
        func_dict['src_var_name'] = var
        init_processing_func(func_dict)
    return var_processing_dict


def process_and_add(gldata, var_name, func_dict):
    """Run any processing functions on the glider data

    :param gldata: GliderData instance to add the processed variable to
    :param var_name: string
        Variable Name, the name of the variable that will be added to the
        GliderData instance which will ultimately become a netCDF variable
    :param func_dict:
    :return: gdata: GliderData instance with new variable added
    """

    # func_dict = init_processing_func(func_dict)
    func = func_dict['func']
    func_name = func.__name__

    # assemble the arguments from the func_dict
    # first get the actual data arguments

    try:
        args = deepcopy(func_dict['data_args'])
    except KeyError:
        logger.error(
            "Processing `{:s}` failed.\n"
            "There MUST be a 'data_args' dictionary in the processing JSON "
            "block that has at least one argument to the function {:s} "
            "whose value is a variable name from the data file.".format(
                var_name, func_name))
        raise

    for arg in args:
        srcvarname = args[arg]
        assert isinstance(srcvarname, str), (
            "Arguments in the 'data_args' dictionary in the {:s} JSON "
            "processing block must be a string variable name from the "
            "glider source data".format(var_name))
        try:
            # return a deepcopy here to prevent a user provided function from
            # accidentally changing the raw data
            var = deepcopy(gldata.getdata(srcvarname))
        except SensorError:
            logger.error(
                "During processing steps for `{:s}`, glider variable `{:s}` "
                "could not be found in the data".format(var_name, srcvarname))
            raise
        else:
            args[arg] = var
    # then add any parameter arguments
    if 'param_args' in func_dict:
        args.update(func_dict['param_args'])

    # noinspection PyBroadException
    # ^ this turns off the Pycharm warning for using Exception as too broad
    # ToDo: these logger error statements should be refactored to just change
    #  the text in the error itself rather than logging, if you want to
    #  capture the error in logging specifically (this program only outputs
    #  the messages to stdout and stderr, so it is a moot point) then you can
    #  use logger.exception instead (the moot point is that it would print
    #  the error message twice to the screen.
    try:
        outvar = func(**args)
    except TypeError:
        logger.error(
            "During processing for `{:s}`, the arguments for function "
            "`{:s}` are incorrect. Check the inputs in the json processing "
            "dictionary. See instructions/how_to_process.txt.".format(
                var_name, func_name))
        raise
    except Exception:
        logger.error(
            "During processing for `{:s}`, the function `{:s}` failed. Test "
            "the function and make sure that it returns an 1D numpy.array the "
            "same length as the input data.".format(var_name, func_name)
        )
        raise

    # ToDo: This section will not work if the function returns multiple
    #  output variables, each the same length as the glider data
    #  For now I want to complete this for processing, but in a future
    #  implementation I want to tackle the challenging scenarios of multiple
    #  outputs and processing from an already processed variable.

    if len(outvar) != len(gldata):
        error_msg = (
            "The length of the output variable from the `{:s}`: `{:s}` "
            "JSON processing block does not match the input data length"
            "Test the processing function's output".format(var_name, func_name)
        )
        logging.error(error_msg)
        raise ValueError(error_msg)

    gldata.add_data(var_name, np.array(outvar))

    # not complete yet.  need to account for the case of isinstance np.array
    # where it is just one variable
    return gldata