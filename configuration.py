

# list the science variables/sensor names from the raw glider data to process.
# These should match with the netCDF variables in sensor_defs.json

# ToDo: incorporate the constants from ooidac/constants.py into here and then
#  setup a system where a required sensor is either a single string or an
#  order of preference list of strings where it will look for the first
#  preference in the data and if not found, go down the list until one is
#  found or if none, then fail or continue.

# ToDo: rename all *SENSOR variables to *VARNAME or *VARIABLE etc. for a more
#  generic approach.  Only Slocum calls variables "sensors" EVERYone else
#  calls them "variables"

# raw glider sensor/variable names
TIMESENSOR = 'm_present_time'
SCITIMESENSOR = 'sci_m_present_time'
CONDUCTIVITY_SENSOR = 'sci_water_cond'
TEMPERATURE_SENSOR = 'sci_water_temp'
PRESSURE_SENSOR = "sci_water_pressure"
DEPTH_SENSOR = 'm_depth'
LAT_SENSOR = "m_gps_lat"
LON_SENSOR = "m_gps_lon"

# processed sensor/variable names
PROC_LAT_VAR = "llat_latitude"
PROC_LON_VAR = "llat_longitude"
PROC_TIME_VAR = "llat_time"
PROC_PRES_VAR = "llat_pressure"
PROC_DEPTH_VAR = "depth"

# netCDF variable names
NC_TIME_VAR = "time"
NC_LON_VAR = "lon"
NC_LAT_VAR = "lat"
NC_PRES_VAR = "pressure"
NC_DEPTH_VAR = "depth"

REQUIRED_SENSORS = [
    SCITIMESENSOR, TIMESENSOR,
    CONDUCTIVITY_SENSOR, TEMPERATURE_SENSOR,
    PRESSURE_SENSOR,
    LAT_SENSOR, LON_SENSOR,
]

# The Depth Averaged Velocity (DAV) sensors/variables
DAV_SENSORS = [
    ('m_final_water_vx', 'm_final_water_vy'),  # first preference
    ('m_initial_water_vx', 'm_initial_water_vy'),  # second preference
    ('m_water_vx', 'm_water_vy')  # third preference
]

# This should be the list of sensors/variables that provide scientific data
DATA_CONFIG_LIST = [
    TEMPERATURE_SENSOR, CONDUCTIVITY_SENSOR, PRESSURE_SENSOR,
    'sci_flbbcd_chlor_units', 'sci_flbbcd_cdom_units', 'sci_flbbcd_bb_units',
    'sci_oxy4_oxygen', 'sci_oxy4_saturation', 'sci_bsipar_par'
]

# The minimum required number of measurements for a variable/sensor to be
# considered valid in a segment data file
MIN_DATA_VALS = 5

# The minimum required depth for the glider to dive to consider a profile or
# a data file to be valid
MIN_DIVE_DEPTH = 2.0


# PROCESSING_DIR = "C:\\Users\\spearce\\data\\dac\\processing_test"
