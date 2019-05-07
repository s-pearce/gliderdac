#!/usr/bin/env python

import logging
import os
import argparse
import sys
import glob

sys.path.append("C:\\Users\\spearce\\code\\python\\gliderdac")

import numpy as np
from netCDF4 import num2date
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import ooidac.processing as process_dba
from ooidac.data_classes import DbaData
from ooidac.profiles import Profiles
from makencw import makencw

_days_formatter = mdates.DayLocator()
_hours_formatter = mdates.HourLocator()
_days_format = mdates.DateFormatter('%m/%d')
_hours_format = mdates.DateFormatter('%H:%M')


def main(args):
    """Parse the specified Slocum glider dba file and plot the yo and the
    indexed profiles"""

    # Set up logger
    log_level = getattr(logging, args.loglevel.upper())
    log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
    logging.basicConfig(format=log_format, level=log_level)

    dba_files = glob.glob(args.dba_files)
    ncw = makencw()  # ToDo: replace the call here with the long term solution
    ctd_sensors = [
        'llat_latitude', 'llat_longitude',
        'llat_pressure', 'sci_water_temp',
        'sci_water_cond']

    for dba_file in dba_files:

        if not os.path.isfile(dba_file):
            logging.error('Invalid dba file specified: {:s}'.format(
                dba_file))
            continue

        logging.debug('Parsing dba file: {:s}'.format(dba_file))
        dba = DbaData(dba_file)
        if dba is None or len(dba) == 0:
            logging.warning('Empty dba file: {:s}'.format(dba_file))
            continue
        dba = process_dba.create_llat_sensors(dba)
        if dba is None:
            continue

        dba = process_dba.ctd_data(dba, ctd_sensors, ncw)
        if dba is None:
            continue

        # Create profile indexing
        logging.debug('Indexing profiles...')
        profiles = Profiles(dba)
        profiles.find_profiles_by_depth(10)
        if len(profiles) == 0:
            logging.info('No profiles indexed: {:s}'.format(dba_file))
            continue

        if args.clean:
            logging.debug('Cleaning up indexed profiles...')
            profiles.filter_profiles()

        if len(profiles) == 0:
            continue

        logging.debug('Plotting profiles...')
        plot_profiles(profiles)
        inflections = num2date(
            profiles.inflection_times,
            'seconds since 1970-01-01 00:00:00Z'
        )
        plot_vert_lines(inflections)
        plt.show()

    return 0


def plot_profiles(profiles):
    """Plot the glider yo and the indexed profiles"""

    ax = plt.gca()

    ax.xaxis.set_major_locator(_days_formatter)
    ax.xaxis.set_major_formatter(_days_format)
    ax.xaxis.set_minor_locator(_hours_formatter)
    ax.xaxis.set_minor_formatter(_hours_format)

    seg_depth = profiles.dba.getdata('m_depth')
    ts = num2date(
        profiles.dba.getdata('llat_time'),
        'seconds since 1970-01-01 00:00:00Z'
    )

    plt.plot(ts, seg_depth, 'k.', axes=ax)
    ax.invert_yaxis()
    colorstyle = [
        'b', 'r', 'g',
        'xkcd:orange', 'xkcd:purple', 'xkcd:yellow',
        'xkcd:pink', 'xkcd:light neon green', 'xkcd:neon red'
    ]
    colorgen = (colorstyle[ii % (len(colorstyle)-1)] for ii in range(len(
        profiles)))

    for profile in profiles:

        p_times = num2date(
            profile.getdata('llat_time'),
            'seconds since 1970-01-01 00:00:00Z'
        )
        p_data = profile.getdata('m_depth')
        color = next(colorgen)
        plt.plot(
            p_times, p_data,
            marker='^', markerfacecolor=color,
            markeredgecolor=color, linestyle=""
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


if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser(
        description=main.__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument(
        'dba_files',
        help='Path to a dba file or files')

    arg_parser.add_argument(
        '-c', '--clean',
        help='Clean up/filter the indexed profiles',
        action='store_true')

    arg_parser.add_argument(
        '-l', '--loglevel',
        help='Verbosity level',
        type=str,
        choices=['debug', 'info', 'warning', 'error', 'critical'],
        default='info')

    parsed_args = arg_parser.parse_args()

    # print(parsed_args)
    # sys.exit(1)

    sys.exit(main(parsed_args))