import os
import glob
import logging
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from configuration import DEPTH_SENSOR, PRESSURE_SENSOR
from ooidac.plotting import color_gen2
from pandas.plotting import register_matplotlib_converters
import netCDF4
from ooidac.data_classes import DbaData
from ooidac.utilities import date_xlabel

register_matplotlib_converters()

REGEX = re.compile(r'[A-Za-z0-9_]+_(\d{8}T\d{6})Z_(?:delayed|rt).nc')


def plot_dba(ts, depth):
    plt.plot(ts, depth, 'k.')
    ax = plt.gca()
    ax.invert_yaxis()
    plt.ylabel("Depth, m")
    plt.xlabel(date_xlabel(ts))
    # plt.xlabel("Time")


def get_profile_graph_angle(profile_time, profile_depth, ax):
    # Get the DATA slope of the line, which is in data units not in on screen
    # figure units.
    # Note: must have no NaN entries
    finites = np.isfinite(profile_depth)
    slope = np.polyfit(profile_time[finites], profile_depth[finites], 1)[0]

    # Convert slope of meters/secs to meters/days since Matplotlib.dates are
    # handled as serial days since 0001-01-01 + 1.  See matplotlib.dates
    # description.
    slope *= 86400  # m/secs * secs/day = m/day

    # Total figure size
    fig_w, fig_h = ax.get_figure().get_size_inches()
    # Axis size on figure
    _, _, w, h = ax.get_position().bounds
    # Ratio of display units
    y_axis_size = (fig_h * h)
    x_axis_size = (fig_w * w)
    disp_ratio = y_axis_size / x_axis_size
    # Ratio of data units
    # Negative over negative because of the order of subtraction
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    dxlim = xlim[1] - xlim[0]
    dylim = ylim[1] - ylim[0]

    data_ratio = dylim / dxlim

    display_line_slope = slope * disp_ratio / data_ratio
    display_line_angle = np.degrees(np.arctan(display_line_slope))

    return display_line_angle


def profile_text_line(ptime, pdepth, msg, ax):
    median_depth = np.nanmedian(pdepth)
    minimizer = abs(pdepth - median_depth)
    pmid_ii = np.flatnonzero(minimizer == np.nanmin(minimizer))[0]
    pmid_depth = pdepth[pmid_ii]
    pmid_time = ptime[pmid_ii]

    ptime_ts = ptime.astype(np.int64) / 1e9
    line_angle = get_profile_graph_angle(ptime_ts, pdepth, ax)

    # I want to offset the text position by 7% of the profile's display
    # width so that it sits next to the line
    p_time_tdelta = ptime[-1] - ptime[0]
    text_x = pmid_time + p_time_tdelta * .07
    plt.text(
        text_x, pmid_depth, msg, rotation=line_angle,
        verticalalignment='center', horizontalalignment='center')


def plot_multiprofiles_and_dba(nc_file_list, plot_dir=None):
    """Takes a list of netCDF files corresponding to the profiles of a glider
    segment plotting the profiles and uses metdata from the NetCDF profiles to
    get the source ascii glider segment data file (must be on the same
    machine) and plotting the original depth vs time plot to compare profile
    discovery to the original file as a visual check for accuracy.

    :param nc_file_list: List of DAC NetCDF Profile files
    :param plot_dir: Directory where plotting are to be saved.  If `None` (
    default) than plotting are printed to the screen one at a time (the next is
    generated upon the closing of a previous figure)
    :return:
    """
    dba_path = None
    dba = None
    source_file = "delete_me"
    title = ""
    styles = color_gen2()
    for ncfile in nc_file_list:
        gldata = netCDF4.Dataset(ncfile)
        prof_dba_path = gldata.variables['source_file'].full_path
        ncbasename = os.path.basename(ncfile)
        match = REGEX.match(ncbasename)
        if match is None:
            logging.error('ncfile {:s} is not named correctly'.format(
                ncbasename))
            return
        ncfile_timestr = match.group(1)
        logging.debug("Plotting profile {:s}".format(ncfile))
        if prof_dba_path != dba_path:
            if dba and plot_dir:
                fig_file = os.path.join(
                    plot_dir, "{:s}.png".format(
                        source_file))
                logging.debug("saving figure as {:s}".format(fig_file))
                plt.savefig(fig_file, dpi=300)
            else:
                plt.show()
            plt.clf()
            logging.info("opening segment {:s}".format(prof_dba_path))
            dba = DbaData(prof_dba_path)
            dba_path = dba.file_metadata['full_path']
            source_file = dba.file_metadata['source_file']
            logging.debug("plotting segment {:s}".format(source_file))
            dba_ts = pd.to_datetime(dba.ts * 1e9)
            # dba_depth = dba.depth
            dba_depth = dba.getdata(DEPTH_SENSOR)
            dba_press = dba.getdata(PRESSURE_SENSOR) * 10.
            plot_dba(dba_ts, dba_depth)
            plt.plot(dba_ts, dba_press, 'gray',
                     marker='.', linestyle="none", markersize=5)
            title = "{:s}:".format(source_file)
            styles = color_gen2()

        prof_ts = pd.to_datetime(gldata['time'][:].data * 1e9)
        prof_depth = gldata['depth'][:].data
        color = next(styles)
        plt.plot(prof_ts, prof_depth, color, marker='.', markersize=2)
        plt.axvspan(prof_ts[0], prof_ts[-1], alpha=0.3, color=color)
        ax = plt.gca()
        profile_text_line(prof_ts, prof_depth, ncfile_timestr, ax)

        plt.title(title)
    if dba and plot_dir:
        fig_file = os.path.join(
            plot_dir, "{:s}.png".format(
                source_file))
        plt.savefig(fig_file, dpi=300)
    else:
        plt.show()


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        description="""Create plotting to verify profiles compared to segment 
        depths""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    arg_parser.add_argument('nc_files',
                            help='Profile netCDF files to check',
                            nargs='+')

    arg_parser.add_argument('-s', '--save_path',
                            help=(
                                'Save figures to the plot directory instead of '
                                'plotting on screen. If no directory is '
                                'given, the plotting will appear one at a time '
                                'on screen.  Closing one plot will then '
                                'generate the next.'))

    arg_parser.add_argument('-l', '--loglevel',
                            help='Verbosity level',
                            type=str,
                            choices=[
                                'debug', 'info', 'warning',
                                'error', 'critical'],
                            default='info')

    parsed_args = arg_parser.parse_args()

    if parsed_args.save_path is not None:
        plot_dir_path = parsed_args.save_path
    else:
        plot_dir_path = None

    if len(parsed_args.nc_files) == 1:
        nc_files = glob.glob(parsed_args.nc_files[0])
    else:
        nc_files = parsed_args.nc_files

    log_level = parsed_args.loglevel

    log_format = '%(levelname)s: [line %(lineno)d] %(message)s'
    log_level = getattr(logging, log_level.upper())
    logging.basicConfig(format=log_format, level=log_level)

    plot_multiprofiles_and_dba(nc_files, plot_dir_path)
