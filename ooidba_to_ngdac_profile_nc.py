#!/usr/bin/env python

import os
import sys
import logging
import argparse
import shutil
# import pdb
import glob
# temporary addition to test this script -SP 2019-01-30 to be able to load
# gncutils -SP 2019-03-11 Maybe don't need after all.  We'll see.
# sys.path.append('C:\\Users\\spearce\\code\\python\\gliderdac\\')
from ooidac.writers.netCDFwriter import NetCDFWriter

from ooidac.constants import NETCDF_FORMATS, LLAT_SENSORS
from ooidac.validate import validate_sensors, validate_ngdac_var_names

import ooidac.processing as processing
from ooidac.data_classes import DbaData
from ooidac.profiles import Profiles
from ooidac.data_checks import check_file_goodness, check_for_dav_sensors
from ooidac.constants import SCI_CTD_SENSORS
from dba_file_sorter import sort_function


def main(args):
    """Parse one or more Slocum glider ascii dba files and write IOOS
    NGDAC-compliant Profile NetCDF files
    """

    # Set up logger, see LogManager class below.  This is used to create
    # different logging message formats, one to announce running the script
    # for a data file, and an indented format to show all log messages that
    # occur under that data file.  This is implemented by
    # LogManager.update_format
    start_log_format = (
        '%(levelname)s:%(module)s: [line %(lineno)d]'
        '\n%(message)s')
    run_log_format = (
        '    %(levelname)s:%(module)s: [line %(lineno)d]'
        '\n        %(message)s')
    logmanager = LogManager(start_log_format, args.loglevel)

    # Gather configurations from the arguments
    config_path = args.config_path
    output_path = args.output_path or os.path.realpath(os.curdir)
    start_profile_id = args.start_profile_id
    clobber = args.clobber
    comp_level = args.compression
    nc_format = args.nc_format
    ctd_sensor_prefix = args.ctd_sensor_prefix

    ctd_sensor_names = [
        '{:s}_{:s}'.format(ctd_sensor_prefix, ctd_sensor) for ctd_sensor in
        ['water_cond', 'water_temp']
    ]
    ctd_sensors = LLAT_SENSORS + ctd_sensor_names
    
    # added the ability to use glob expansion; SP Feb2019
    if len(args.dba_files) == 1:
        dba_files = glob.glob(args.dba_files[0])
    else:
        dba_files = args.dba_files
    # sort dba files since the dbd format doesn't sort well.
    dba_files.sort(key=sort_function)

    # Initial Checks
    if not os.path.isdir(config_path):
        logging.error(
            'Invalid configuration directory: {:s}'.format(config_path))
        return 1

    # TODO: checked.  this is mildly redundant with the lines above.  So
    #  remove redundancy and fix.
    if not output_path:
        args.output_path = os.path.realpath(os.curdir)
        logging.info(
            'No NetCDF output_path specified. Using cwd: {:s}'.format(
                output_path))

    if not os.path.isdir(output_path):
        logging.error('Invalid output_path: {:s}'.format(output_path))
        return 1

    if not dba_files:
        logging.error('No Slocum dba files specified')
        return 1

    # Create the Trajectory NetCDF writer
    ncw = NetCDFWriter(
        config_path, output_path, comp_level=comp_level,
        nc_format=nc_format, profile_id=start_profile_id,
        clobber=clobber)
    # Make sure we have llat_* sensors defined in ncw.nc_sensor_defs
    ctd_valid = validate_sensors(ncw.nc_sensor_defs, ctd_sensors)
    if not ctd_valid:
        logging.error(
            'Bad sensor definitions: {:s}'.format(ncw.sensor_defs_file))
        return 1
    # Make sure we have configured sensor definitions for all IOOS NGDAC
    # required variables
    ngdac_valid = validate_ngdac_var_names(ncw.nc_sensor_defs)
    if not ngdac_valid:
        logging.error(
            'Bad sensor definitions: {:s}'.format(ncw.sensor_defs_file))
        return 1

    if args.debug:
        sys.stdout.write('{}\n'.format(ncw))
        return 0

    # Create a temporary directory for creating/writing NetCDF prior to
    # moving them to output_path

    # update the logging format so that indentation can show log statements
    # sub-level to the file being processed after an initial processing file
    # statement

    # Write one NetCDF file for each input file
    output_nc_files = []
    processed_dbas = []
    for dba_file in dba_files:
        # change to non-indented log format (see above)
        logmanager.update_format(start_log_format)

        if not os.path.isfile(dba_file):
            logging.error('Invalid dba file specified: {:s}'.format(dba_file))
            continue

        logging.info('Processing data file: {:s}'.format(dba_file))

        # change to indented log format (see above)
        logmanager.update_format(run_log_format)

        # Parse the dba file
        dba = DbaData(dba_file)

        if dba is None or dba.N == 0:
            logging.warning('Skipping empty data file: {:s}'.format(dba_file))
            continue
        mission = dba.file_metadata['mission_name'].upper()
        if (
                mission == 'STATUS.MI'
                or mission == 'LASTGASP.MI'
                or mission == 'INITIAL.MI'):
            logging.info('Skipping {:s} data file'.format(mission))
            continue

        # check the data file for the required sensors and that science data
        # exists in the file.  True or False returned.  See data_checks.py
        file_good = check_file_goodness(dba)
        if not file_good or len(dba.underwater_indices) == 0:
            logging.warning(
                'File {:s} either does not have enough science data or lacks '
                'the required sensors to produce DAC formatted profiles'.format(
                    dba.source_file)
            )
            continue

        scalars = []

        # remove any sci bay intialization zeros that may occur (where all
        # science instrument sensors are 0.0)
        dba = processing.remove_initial_sci_zeros(dba)

        # This processing step adds time dependent coordinate variables
        # (designated by the prefix llat [lat, lon, altitude, time]) which are
        # created and added to the data instance with metadata attributes.
        # This step is required before the other processing module steps.  The
        # variable llat_time is derived from `m_present_time`,
        # llat_latitude/longitude are filled by linear interpolation from
        # `m_gps_lat/lon`, llat_pressure is derived from converting
        # `sci_water_pressure` to dbar, and llat_depth is derived from
        # `llat_pressure` converted to depth using the Python TEOS-10 GSW
        # package
        dba = processing.create_llat_sensors(dba)
        if dba is None:
            continue

        # Convert m_pitch and m_roll variables to degrees, and add back to
        # the data instance with metadata attributes
        if 'm_pitch' in dba.sensor_names and 'm_roll' in dba.sensor_names:
            dba = processing.pitch_and_roll(dba)
            if dba is None:
                continue

        # Convert `sci_water_cond/temp/ & pressure` to `salinity` and `density`
        # and adds them back to the data instance with metadata attributes.
        # Requires `llat_latitude/longitude` variables are in the data
        # instance from the `create_llat_sensors` method
        dba = processing.ctd_data(dba, SCI_CTD_SENSORS)
        if dba is None:
            continue

        # Process `sci_oxy4_oxygen` to OOI L2 compensated for salinity and
        # pressure and converted to umol/kg.
        if 'sci_oxy4_oxygen' in dba.sensor_names:
            dba = processing.o2_s_and_p_comp(dba)
            if dba is None:
                continue

        if 'sci_flbbcd_bb_units' in dba.sensor_names:
            dba = processing.backscatter_total(dba)
            if dba is None:
                continue
            radiation_wavelength = {
                'data': 700,
                'attrs': {'units': 'nm'},
                'nc_var_name': 'radiation_wavelength'}
            scalars.append(radiation_wavelength)

        # If Depth Averaged Velocity (DAV) data available, (i.e. any of the
        # `*_water_vx/vy` sensors are in the data) get the values and calculate
        # the mean position and time of the segment as the postion and time
        # for the velocity estimation
        if check_for_dav_sensors(dba) and len(dba.underwater_indices) > 0:
            # get segment mean time, segment mean lat, and segment mean lon
            # (underwater portion only)
            seg_time, seg_lat, seg_lon = processing.get_segment_time_and_pos(
                dba)
            if seg_time is None or seg_lat is None or seg_lon is None:
                break

            # if data is the recovered data, this tries to get
            # `m_final_water_vx/vy` from the next 2 segment data files where
            # the calculation occurs, if it cannot, it will get either
            # `m_initial_water_vx/vy` or `m_water_vx/vy` from the current
            # segement file.
            dba_index = dba_files.index(dba_file)
            next2files = dba_files[dba_index + 1:dba_index + 3]
            vx, vy = processing.get_u_and_v(dba, check_files=next2files)

            scalars.extend([seg_time, seg_lat, seg_lon, vx, vy])

        # create a config object
        # create a nc_defs object
        # create a global attributes object
        # create an instruments object
        # create a deployment object

        # The Profiles class discovers the profiles and filters them based on
        # user configurable filters written in profile_filters.py (Instructions
        # found in the module) and returns each profile as a new GliderData
        # instance that is the data subset of the main data instance
        profiles = Profiles(dba)

        profiles.find_profiles_by_depth()

        # See profile_filters.py for which filters are applied
        profiles.filter_profiles()

        if len(profiles.indices) == 0:
            logging.info('No profiles indexed: {:s}'.format(dba_file))
            continue

        # TODO: can move this check to where lats and lons are processed
        #  Note: this is actually done in process_dba.process_ctd, but maybe
        #  it should also be done where lat and lon is processed
        # Make sure we have latitudes and longitudes
        # if np.all(np.isnan(ctd_data[:, 2])):
        #     logging.warning(
        #         'dba contains no valid llat_latitude values'.format(dba_file))
        #     logging.info('Skipping dba: {:s}'.format(dba_file))
        #     continue
        # if np.all(np.isnan(ctd_data[:, 3])):
        #     logging.warning(
        #         'dba contains no valid llat_longitude values'.format(
        #         dba_file))
        #     logging.info('Skipping dba: {:s}'.format(dba_file))
        #     continue

        # %------ Write Profiles to NetCDF ------%
        for profile in profiles:
            profile = processing.reduce_to_sci_data(profile)
            # ToDo: fix the history writer in NetCDFWriter
            out_nc_file = ncw.write_profile(profile, scalars)
            if out_nc_file:  # can be None if skipping
                output_nc_files.append((dba_file, out_nc_file))

        processed_dbas.append(dba_file)

    # Delete the temporary directory once files have been moved
    try:
        logging.debug('Removing temporary files:')
        shutil.rmtree(ncw.tmp_dir)
    except OSError as e:
        logging.error(e)
        logging.error('manually try to remove temporary directory {:s}'.format(
            ncw.tmp_dir)
        )

    # Print the list of files created
    sys.stdout.write('Profiles NC files written:\n')
    for output_nc_file in output_nc_files:
        base_nc = os.path.basename(output_nc_file[1])
        base_data = os.path.basename(output_nc_file[0])
        os.chmod(output_nc_file[1], 0o664)
        sys.stdout.write('\t{:s} -> {:s}\n'.format(base_data, base_nc))

    return 0


class LogManager:

    def __init__(self, log_format, log_level, **kwargs):
        self.log_level = getattr(logging, log_level.upper())
        logging.basicConfig(level=self.log_level, **kwargs)
        self.logger = logging.getLogger()
        self.logger_handler = self.logger.handlers[0]
        self.update_format(log_format)

    def update_format(self, msg_format):
        self.logger_handler.setFormatter(logging.Formatter(msg_format))

    def get_logger(self):
        return self.logger


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(
        description=main.__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    arg_parser.add_argument('config_path',
                            help='Location of deployment configuration files')

    arg_parser.add_argument('dba_files',
                            help='Source ASCII dba files to process',
                            nargs='+')

    arg_parser.add_argument('--ctd_sensor_prefix',
                            help=(
                                'Specify Slocum glider ctd sensor prefix '
                                'letter, i.e.: sci_water_temp, sci_water_cond'),
                            choices=['sci', 'm'],
                            default='sci')

    arg_parser.add_argument('-p', '--start_profile_id',
                            help=(
                                'Integer specifying the beginning profile '
                                'id. If not specified or <1 the mean profile '
                                'unix timestamp is used'),
                            type=int,
                            default=0)

    arg_parser.add_argument('-o', '--output_path',
                            help=(
                                'NetCDF destination directory, which must '
                                'exist. Current directory if not specified'))

    arg_parser.add_argument('-c', '--clobber',
                            help='Clobber existing NetCDF files if they exist',
                            action='store_true')

    arg_parser.add_argument('-f', '--format',
                            dest='nc_format',
                            help='NetCDF file format',
                            choices=NETCDF_FORMATS,
                            default='NETCDF4_CLASSIC')

    arg_parser.add_argument('--compression',
                            help='NetCDF4 compression level',
                            choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                            default=1)

    arg_parser.add_argument('-x', '--debug',
                            help=(
                                'Check configuration and create NetCDF file '
                                'writer, but does not process any files'),
                            action='store_true')

    arg_parser.add_argument('-l', '--loglevel',
                            help='Verbosity level',
                            type=str,
                            choices=[
                                'debug', 'info', 'warning',
                                'error', 'critical'],
                            default='info')

    parsed_args = arg_parser.parse_args()

    sys.exit(main(parsed_args))
