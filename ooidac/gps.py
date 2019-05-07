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


def get_decimal_degrees(lat_lon):
    """Converts glider gps coordinate ddmm.mmm to decimal degrees dd.ddd

    Arguments:
    lat_lon - A floating point latitude or longitude in the format ddmm.mmm
        where dd's are degrees and mm.mmm is decimal minutes.

    Returns decimal degrees float
    """

    if np.isnan(lat_lon):
        return np.nan
        
    # Absolute value of the coordinate
    try:
        pos_lat_lon = abs(lat_lon)
    except (TypeError, ValueError) as e:
        return
    
    # Calculate NMEA degrees as an integer
    nmea_degrees = int(pos_lat_lon/100)*100
    
    # Subtract the NMEA degrees from the absolute value of lat_lon and divide
    # by 60 to get the minutes in decimal format
    gps_decimal_minutes = (pos_lat_lon - nmea_degrees)/60.0
    
    # Divide NMEA degrees by 100 and add the decimal minutes
    decimal_degrees = (nmea_degrees/100) + gps_decimal_minutes
    
    if lat_lon < 0:
        return -decimal_degrees
    
    return decimal_degrees


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
