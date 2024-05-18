#!/usr/bin/env python
import logging

import numpy as np
from gsw import SP_from_C, SA_from_SP, CT_from_t, rho, z_from_p
from gsw.density import sound_speed
from ooidac.constants import SLOCUM_REALTIME_MODE_EXTENSIONS
from ooidac.processing import all_sci_indices

def create_practical_salinity_sensor(reader):

    return


def create_density_sensor(reader):

    return


def calculate_practical_salinity(conductivity, temperature, pressure):
    """Calculates practical salinity given glider conductivity, temperature,
    and pressure using Gibbs gsw SP_from_C function.

    Parameters:
        timestamp, conductivity (S/m), temperature (C), and pressure (bar).

    Returns:
        salinity (psu PSS-78).
    """

    # Convert S/m to mS/cm
    ms_conductivity = conductivity * 10

    return SP_from_C(
        ms_conductivity,
        temperature,
        pressure
    )


def calculate_density(temperature, pressure, salinity, latitude, longitude):
    """Calculates density given glider practical salinity, pressure, latitude,
    and longitude using Gibbs gsw SA_from_SP and rho functions.

    Parameters:
        timestamps (UNIX epoch),
        temperature (C), pressure (dbar), salinity (psu PSS-78),
        latitude (decimal degrees), longitude (decimal degrees)

    Returns:
        density (kg/m**3),
    """

    # dBar_pressure = pressure * 10

    absolute_salinity = SA_from_SP(
        salinity,
        pressure,
        longitude,
        latitude
    )

    conservative_temperature = CT_from_t(
        absolute_salinity,
        temperature,
        pressure
    )

    density = rho(
        absolute_salinity,
        conservative_temperature,
        pressure
    )

    return density


def calculate_sound_speed(temperature, pressure, salinity, latitude, longitude):
    """Calculates sound speed given glider practical in-situ temperature, pressure and salinity using Gibbs gsw
    SA_from_SP and rho functions.

    Parameters:
        temperature (C), pressure (dbar), salinity (psu PSS-78), latitude, longitude

    Returns:
        sound speed (m s-1)
    """

    absolute_salinity = SA_from_SP(
        salinity,
        pressure,
        longitude,
        latitude
    )

    conservative_temperature = CT_from_t(
        absolute_salinity,
        temperature,
        pressure
    )

    speed = sound_speed(
        absolute_salinity,
        conservative_temperature,
        pressure
    )

    return speed


def calculate_depth(pressure, latitude):
    """Calculates depth from pressure (dbar) and latitude.  By default, gsw returns depths as negative.  This routine
    returns the absolute values for positive depths.

    Parameters:
        pressure (decibars)
        latitude (decimal degrees)

    Returns:
        depth (meters)
    """

    return abs(z_from_p(pressure, latitude))


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

    # There is a bug in the tbdlist.dat decimation process on the gliders
    # that occasionally misaligns the CTD variables in time on TBD files.
    # This isn't a problem in EBD files.  This if statement interpolates
    # to correct this.  For ease of use, it will just overwrite the original
    # 'sci_water_temp', 'sci_water_cond', and 'sci_ctd41cp_timestamp'.  No
    # need to overwrite 'sci_water_pressure' because we use the coordinate
    # pressure which is already interpolated to all timestamps.
    ext = dba.file_metadata['filename_extension']
    if ext in SLOCUM_REALTIME_MODE_EXTENSIONS:
        # I'm just going to go ahead and assume if ctd data is present,
        # these variables are there too.  I know this is bad practice,
        # but I plan to change most of this eventually anyway.
        rawpres = dba['sci_water_pressure'] or dba['m_water_pressure']
        sts = dba['sci_m_present_time']['data']

        cond_ii = np.isfinite(cond['data'])
        temp_ii = np.isfinite(temp['data'])
        pres_ii = np.isfinite(rawpres['data'])

        ctd_indices = np.logical_or(cond_ii, np.logical_or(temp_ii, pres_ii))

        icond = np.full_like(cond['data'], np.nan)
        itemp = np.full_like(temp['data'], np.nan)

        icond[ctd_indices] = np.interp(
            sts[ctd_indices], sts[cond_ii], cond['data'][cond_ii])
        itemp[ctd_indices] = np.interp(
            sts[ctd_indices], sts[temp_ii], temp['data'][temp_ii])

        dba['sci_water_cond']['data'] = icond
        dba['sci_water_cond']['attrs']['comment'] = (
            "Interpolated to all CTD timestamps to re-align "
            "realtime file decimation misalignments")

        dba['sci_water_temp']['data'] = itemp
        dba['sci_water_temp']['attrs']['comment'] = (
            "Interpolated to all CTD timestamps to re-align "
            "realtime file decimation misalignments")

        # if the ctd timestamp is present, interpolate it too
        if 'sci_ctd41cp_timestamp' in dba.sensor_names:
            tctd = dba['sci_ctd41cp_timestamp']
            tctd_ii = np.isfinite(tctd['data'])
            itctd = np.full_like(tctd['data'], np.nan)
            itctd[ctd_indices] = np.interp(
                sts[ctd_indices], sts[tctd_ii], tctd['data'][tctd_ii])
            tctd['data'] = itctd

    # Calculate mean llat_latitude and mean llat_longitude
    mean_lat = np.nanmean(lat['data'])
    mean_lon = np.nanmean(lon['data'])

    # Calculate practical salinity
    prac_sal = calculate_practical_salinity(
        cond['data'], temp['data'], pres['data'])
    # Add salinity to the dba
    dba['salinity'] = {
        'sensor_name': 'salinity',
        'attrs': {
            "source_sensor": "sci_water_cond,sci_water_temp,pressure"
        },  # these get filled in later by the netCDF writer
        'data': prac_sal
    }

    # Calculate density
    density = calculate_density(
        temp['data'], pres['data'], prac_sal, mean_lat, mean_lon)
    # Add density to the dba
    dba['density'] = {
        'sensor_name': 'density',
        'attrs': {
            "source_sensor": "sci_water_cond,sci_water_temp,pressure,latitude,longitude"
        },  # these get filled in later by the netCDF writer
        'data': density
    }
    return dba


class CTDprocessingError(Exception):
    pass