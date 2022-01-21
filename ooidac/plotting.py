import os
import logging
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from ooidac.processing.ctd import ctd_data

from configuration import PROC_LON_VAR, PROC_LAT_VAR, PROC_PRES_VAR, PROC_TIME_VAR
from configuration import CONDUCTIVITY_SENSOR, TEMPERATURE_SENSOR, DEPTH_SENSOR
from configuration import PRESSURE_SENSOR
import ooidac.processing as process_dba
from ooidac.data_classes import DbaData
from ooidac.profiles import Profiles

_days_formatter = mdates.DayLocator()
_hours_formatter = mdates.HourLocator()
_days_format = mdates.DateFormatter('%m/%d')
_hours_format = mdates.DateFormatter('%H:%M')

ctd_sensors = [
    PROC_LAT_VAR, PROC_LON_VAR,
    PROC_PRES_VAR, TEMPERATURE_SENSOR,
    CONDUCTIVITY_SENSOR]


def plot_profiles_for_dba(dba_file, unfiltered=False):
    if not os.path.isfile(dba_file):
        logging.error('Invalid dba file specified: {:s}'.format(
            dba_file))
        return

    logging.debug('Parsing dba file: {:s}'.format(dba_file))
    dba = DbaData(dba_file)
    if dba is None or len(dba) == 0:
        logging.warning('Empty dba file: {:s}'.format(dba_file))
        return
    dba = process_dba.create_llat_sensors(dba)
    if dba is None:
        return

    dba = ctd_data(dba, ctd_sensors)
    if dba is None:
        return

    # Create profile indexing
    logging.debug('Indexing profiles...')
    profiles = Profiles(dba)
    profiles.find_profiles_by_depth()

    if len(profiles) == 0:
        logging.info('No profiles indexed: {:s}'.format(dba_file))
        return

    if not unfiltered:
        logging.debug('Filtering indexed profiles...')
        profiles.filter_profiles()

    if len(profiles) == 0:
        logging.warning(
            'All profiles removed by filtering: {:s}'.format(dba_file))
        return

    logging.debug('Plotting profiles...')
    plot_profiles(profiles)
    inflections = list(map(
        dt.datetime.utcfromtimestamp, profiles.inflection_times))

    plot_vert_lines(inflections)
    plt.show()


def plot_profiles(profiles):
    """Plot the glider yo and the indexed profiles"""

    ax = plt.gca()

    ax.xaxis.set_major_locator(_days_formatter)
    ax.xaxis.set_major_formatter(_days_format)
    ax.xaxis.set_minor_locator(_hours_formatter)
    ax.xaxis.set_minor_formatter(_hours_format)

    color_styles = color_gen2()

    seg_depth = profiles.dba.getdata(DEPTH_SENSOR)
    seg_pressure = profiles.dba.getdata(PRESSURE_SENSOR) * 10
    ts = list(map(dt.datetime.utcfromtimestamp, profiles.dba.getdata(
        PROC_TIME_VAR)))

    plt.plot(ts, seg_depth, 'k.', axes=ax)
    plt.plot(ts, seg_pressure, marker='.', color='gray', linestyle="")
    ax.invert_yaxis()

    for profile in profiles:
        p_times = list(map(dt.datetime.utcfromtimestamp, profile.getdata(
            PROC_TIME_VAR)))
        p_data = profile.getdata(PRESSURE_SENSOR) * 10
        color = next(color_styles)
        plt.axvspan(p_times[0], p_times[-1], alpha=0.3, color=color)
        plt.plot(
            p_times, p_data,
            marker='^', linestyle="",
        )
    plt.title(profiles.dba.source_file)


def plot_vert_lines(x, style='k:', ax=None):
    x = np.atleast_1d(x)
    n = len(x)
    x = np.tile(x, (2, 1))
    if not ax:
        ax = plt.gca()
    ylims = ax.get_ylim()
    y = np.array([np.repeat(ylims[0], n), np.repeat(ylims[1], n)])
    plt.plot(x, y, style)


def color_gen():
    color_order = [
        '#9e0168', '#c0fb2d',
        'red', 'green',
        '#fd3c06', '#137e6d',
        'orange', 'blue',
        '#fcb001', '#5d06e9',
        'yellow', 'violet']
    ii = 0
    while True:
        yield color_order[ii % len(color_order)]
        ii += 1


def color_gen2():
    color_order = ['red', 'blue']
    ii = 0
    while True:
        yield color_order[ii % len(color_order)]
        ii += 1