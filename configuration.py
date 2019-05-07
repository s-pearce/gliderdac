

# list the science variables/sensor names from the raw glider data to process.
# These should match with the netCDF variables in sensor_defs.json

conductivity_sensor = 'sci_water_cond'
temperature_sensor = 'sci_water_temp'

REQUIRED_SENSORS = [
    'sci_m_present_time', 'm_present_time',
    conductivity_sensor, temperature_sensor,
    'sci_water_pressure',
    'm_gps_lat', 'm_gps_lon',
]

DATA_CONFIG_LIST = [
    temperature_sensor, conductivity_sensor, 'sci_water_pressure',
    'sci_flbbcd_chlor_units', 'sci_flbbcd_cdom_units', 'sci_flbbcd_bb_units',
    'sci_oxy4_oxygen', 'sci_oxy4_saturation', 'sci_bsipar_par'
]
