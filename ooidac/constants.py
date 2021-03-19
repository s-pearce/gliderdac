import netCDF4
from configuration import PROC_TIME_VAR, PROC_DEPTH_VAR, PROC_PRES_VAR
from configuration import PROC_LON_VAR, PROC_LAT_VAR
from configuration import CONDUCTIVITY_SENSOR, TEMPERATURE_SENSOR, PRESSURE_SENSOR

NC_FILL_VALUES = netCDF4.default_fillvals

# Valid python netCDF4 NetCDF file types
NETCDF_FORMATS = ['NETCDF3_CLASSIC',
                  'NETCDF4_CLASSIC',
                  'NETCDF4']

# Slocum delayed mode file types
SLOCUM_DELAYED_MODE_EXTENSIONS = ['dbd',
                                  'ebd']

# Slocum real-time mode file types
SLOCUM_REALTIME_MODE_EXTENSIONS = ['sbd',
                                   'mbd',
                                   'nbd',
                                   'tbd']

# Slocum native timestamp sensors   
SLOCUM_TIMESTAMP_SENSORS = ['m_present_time',
                            'sci_m_present_time']

# Slocum native pressure sensors
SLOCUM_PRESSURE_SENSORS = ['sci_water_pressure',
                           'm_water_pressure',
                           'm_pressure']

# Slocum native depth sensors
SLOCUM_DEPTH_SENSORS = ['m_depth']

SLOCUM_TEMPERATURE_SENSORS = ['sci_water_temp',
                              'm_water_temp']

SLOCUM_SALINITY_SENSORS = ['sci_water_cond',
                           'm_water_cond']

LLAT_SENSORS = [
    PROC_TIME_VAR, PROC_PRES_VAR, PROC_LAT_VAR, PROC_LON_VAR, PROC_DEPTH_VAR]

NGDAC_VAR_NAMES = ['time', 'depth', 'pressure', 'temperature', 'conductivity',
                   'salinity', 'density', 'profile_time', 'profile_id',
                   'profile_lat', 'profile_lon', 'lat', 'lon']

REQUIRED_SENSOR_DEFS_KEYS = [
    # old: 'nc_var_name',
    'type', 'dimension', 'attrs']
CF_VARIABLE_ATTRIBUTES = ['long_name', 'standard_name', 'units', 'comment']

SCI_CTD_SENSORS = [PROC_LAT_VAR,
                   PROC_LON_VAR,
                   PROC_PRES_VAR,
                   'sci_water_temp',
                   'sci_water_cond']

M_CTD_SENSORS = [PROC_LAT_VAR,
                 PROC_LON_VAR,
                 PROC_PRES_VAR,
                 'm_water_temp',
                 'm_water_cond']
