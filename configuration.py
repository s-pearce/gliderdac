

# list the science variables/sensor names from the raw glider data to process.
# These should match with the netCDF variables in sensor_defs.json

TIMESENSOR = 'm_present_time'
SCITIMESENSOR = 'sci_m_present_time'
CONDUCTIVITY_SENSOR = 'sci_water_cond'
TEMPERATURE_SENSOR = 'sci_water_temp'

REQUIRED_SENSORS = [
    SCITIMESENSOR, TIMESENSOR,
    CONDUCTIVITY_SENSOR, TEMPERATURE_SENSOR,
    'sci_water_pressure',
    'm_gps_lat', 'm_gps_lon',
]

DATA_CONFIG_LIST = [
    TEMPERATURE_SENSOR, CONDUCTIVITY_SENSOR, 'sci_water_pressure',
    'sci_flbbcd_chlor_units', 'sci_flbbcd_cdom_units', 'sci_flbbcd_bb_units',
    'sci_oxy4_oxygen', 'sci_oxy4_saturation', 'sci_bsipar_par'
]

# The minimum required number of measurements for a variable/sensor to be
# considered valid in a segment data file
MIN_DATA_VALS = 5

# The minimum required depth for the glider to dive to consider a profile or
# a data file to be valid
MIN_DIVE_DEPTH = 2.0
