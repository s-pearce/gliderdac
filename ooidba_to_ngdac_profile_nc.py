#!/usr/bin/env python

import os
import sys
import logging
import argparse
import shutil
# import pdb
import glob

import ooidac.processing.attitude
import ooidac.processing.ctd
import ooidac.processing.fluorometer
import ooidac.processing.oxygen
import ooidac.processing.par
import ooidac.processing.velocity
from ooidac.writers.netCDFwriter import NetCDFWriter

from ooidac.constants import NETCDF_FORMATS, LLAT_SENSORS
from ooidac.validate import validate_sensors, validate_ngdac_var_names

import ooidac.processing as processing
from ooidac.data_classes import DbaData
from ooidac.profiles import Profiles
from ooidac.data_checks import check_file_goodness
from ooidac.constants import SCI_CTD_SENSORS
from ooidac.status import Status
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
    else:
        output_path = os.path.realpath(output_path)

    if not os.path.isdir(output_path):
        logging.error('Invalid output_path: {:s}'.format(output_path))
        return 1

    if not dba_files:
        logging.error('No Slocum dba files specified')
        return 1

    # Create a status.json file in the config directory given to hold
    # information from the latest run
    status_path = os.path.join(config_path, 'status.json')

    # if not os.path.exists(status_path):
    #     status = {
    #         "history": "", "date_created": "", "date_modified": "",
    #         "date_issued": "", "version": "", "uuid": "",
    #         "raw_directory": os.path.dirname(os.path.realpath(dba_files[0])),
    #         "nc_directory": output_path,
    #         "next_profile_id": None, "files_processed": [],
    #         "profiles_created": [], "profiles_uploaded": [],
    #         "profile_to_data_map": []
    #         }
    # else:
    #     with open(status_path, 'r') as fid:
    #         status = json.load(fid)
    status = Status(status_path)

    # eliminate files that have already run
    if clobber:
        files_to_run = dba_files
        status.update_modified_date()
        already_processed = []  # see else for description
    else:
        dba_files2 = list(map(os.path.abspath, dba_files))
        files_to_run = set(dba_files2).difference(
            status.info['files_processed'])
        files_to_run = list(files_to_run)
        # save a list of profiles created basenames for printout at the end
        already_processed = list(map(
            os.path.basename, status.info['profiles_created']))

    if files_to_run:
        n_skipped = len(files_to_run) - len(dba_files)
        logging.debug("Skipping {:d} files already run".format(n_skipped))
        files_to_run.sort(key=sort_function)
    else:
        logging.info("No new files to run, exiting.")
        return 0


    # get the next profile id if this dataset has been run before.
    if status.info['next_profile_id'] and start_profile_id > 0 and not clobber:
        start_profile_id = status.info['next_profile_id']


    # Create the Trajectory NetCDF writer
    ncw = NetCDFWriter(
        config_path, output_path, comp_level=comp_level,
        nc_format=nc_format, starting_profile_id=start_profile_id,
        clobber=clobber)
    if not status.trajectory:
        status.trajectory = ncw.trajectory

    # use running status object to write global attributes
    ncw.add_global_attribute(
        'history', status.info['history'], override=True)
    ncw.add_global_attribute(
        'date_created', status.info['date_created'], override=True)
    ncw.add_global_attribute(
        'date_modified', status.info['date_modified'], override=True)
    ncw.add_global_attribute(
        'date_issued', status.info['date_issued'], override=True)
    ncw.add_global_attribute(
        'uuid', status.info['uuid'], override=True)
    ncw.add_global_attribute(
        'version', status.info['version'], override=True)

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

    # update the logging format so that indentation can show log statements
    # sub-level to the file being processed after an initial processing file
    # statement

    # These are now handled by the Status class
    # # Write one NetCDF file for each input file
    # output_nc_files = []
    # source_dba_files = []
    # processed_dbas = []
    # profile_to_data_map = []

    # Pre-processing
    # ToDo: clean this up
    var_processing = {}
    for var_defs in ncw.config_sensor_defs:
        if "processing" in ncw.config_sensor_defs[var_defs]:
            var_processing[var_defs] = ncw.config_sensor_defs[var_defs].pop(
                "processing")
    # to be added:
    # var_processing = processing.init_processing_dict(var_processing)

    for dba_file in files_to_run:
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
            # skip source file in future runs by adding to status
            status.add_src(dba_file)
            continue

        # check the data file for the required sensors and that science data
        # exists in the file.  True or False returned.  See data_checks.py
        file_check = check_file_goodness(dba)
        if not file_check.file_good or len(dba.underwater_indices) == 0:
            logging.warning(
                'File {:s} either does not have enough science data, lacks '
                'the required sensors, or does not have any dives deep enough '
                'to produce DAC formatted profiles'.format(dba.source_file)
            )
            status.add_src(dba_file)
            continue

        scalars = []

        # remove any sci bay intialization zeros that may occur (where all
        # science instrument sensors are 0.0)
        dba = processing.remove_sci_init_zeros(
            dba, file_check.avail_sci_data)

        # for the rare instance where an instrument has been removed from
        # proglets.dat, the variables/sensors associated are removed from the
        # segment data file.  If it is missing and is in the DATA_CONFIG_LIST in
        # the gdac configuration file, then we add it back in here as an
        # array of NaNs.
        dba = processing.replace_missing_sensors(
            dba, file_check.avail_sci_data)

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
            dba = ooidac.processing.attitude.pitch_and_roll(dba)
            if dba is None:
                continue

        # Convert `sci_water_cond/temp/ & pressure` to `salinity` and `density`
        # and adds them back to the data instance with metadata attributes.
        # Requires `llat_latitude/longitude` variables are in the data
        # instance from the `create_llat_sensors` method
        dba = ooidac.processing.ctd.ctd_data(dba, SCI_CTD_SENSORS)
        if dba is None:
            continue

        # Process `sci_oxy4_oxygen` to OOI L2 compensated for salinity and
        # pressure and converted to umol/kg.
        if 'corrected_oxygen' in ncw.config_sensor_defs:
            dba = ooidac.processing.oxygen.check_and_recalc_o2(
                dba,
                calc_type=var_processing['corrected_oxygen'][
                    'calculation_type'],
                cal_dict=var_processing['corrected_oxygen']['cal_coefs']
            )
            dba = ooidac.processing.oxygen.o2_s_and_p_comp(dba, 'temp_corrected_oxygen')
            oxy = dba['oxygen']
            oxy['sensor_name'] = 'corrected_oxygen'
            dba['corrected_oxygen'] = oxy
        elif 'sci_oxy4_oxygen' in dba.sensor_names:
            dba = ooidac.processing.oxygen.o2_s_and_p_comp(dba)
            if dba is None:
                continue

        # Re_calculate chlorophyll
        if 'corrected_chlor' in ncw.config_sensor_defs:
            dba = ooidac.processing.fluorometer.recalc_chlor(
                dba, **var_processing['corrected_chlor']
            )
            if dba is None:
                continue

        # Re_calculate CDOM
        if 'corrected_cdom' in ncw.config_sensor_defs:
            dba = ooidac.processing.fluorometer.recalc_cdom(
                dba, **var_processing['corrected_cdom']
            )
            if dba is None:
                continue

        # Re_calculate PAR
        if 'corrected_par' in ncw.config_sensor_defs:
            # par_sensor_dark = corrections['corrected_par']['sensor_dark']
            # par_sf = corrections['corrected_par']['scale_factor']
            dba = ooidac.processing.par.recalc_par(
                dba, **var_processing['corrected_par']
                # sensor_dark=par_sensor_dark,
                # scale_factor=par_sf
            )
            if dba is None:
                continue

        bksctr_vars = [x for x in ncw.config_sensor_defs if 'backscatter' in x]
        for var_name in bksctr_vars:
            if var_name in var_processing:
                bksctr_args = var_processing[var_name]
                assert "source_sensor" in bksctr_args, (
                    'In the "processing" dictionary for variable "{:s}" in '
                    'sensor_defs.json, a "source_sensor" field must be present '
                    'with the glider sensor name to use for '
                    'processing.'.format(var_name)
                )
                bb_sensor = bksctr_args.get('source_sensor')

            else:
                # for now assume we are using flbbcds with 700 nm wavelength as
                # the default, the input values to the backscatter_total
                # function also default to the flbbcd values.
                bksctr_args = {}
                wavelength = 700.0
                bb_sensor = 'sci_flbbcd_bb_units'

            dba = ooidac.processing.fluorometer.backscatter_total(
                dba, bb_sensor, var_name, **bksctr_args)
            if dba is None:
                continue

        # Add radiation wavelength variables as a paired variable to
        # the backscatter variables
        # if wavelength is not present, assume FLBB 700 nm
        # Note: this will likely change in the future to only include
        #   'radiation_wavelength' as an attribute to the backscatter variable.
        #   It is too much to have a variable for each if there are multiple
        #   backscatter variables.  Although the future release will control if
        #   you want it to be a variable or not just by the processing
        #   dictionary
        rw_vars = [x for x in ncw.config_sensor_defs
                   if 'radiation_wavelength' in x]
        for rw_var in rw_vars:
            sdef = ncw.config_sensor_defs[rw_var]
            if 'radiation_wavelength' in sdef['attrs']:
                wl = sdef['attrs']['radiation_wavelength']
            else:
                wl = 700.0
            radiation_wavelength = {
                'data': wl,
                'attrs': {'units': 'nm'},
                'nc_var_name': rw_var}
            scalars.append(radiation_wavelength)

        # To be added:
        # ---------General Processing------------#
        # if var_processing:
        #     for var in var_processing:
        #         dba = processing.process_and_add(dba, var, var_processing[var])
        #         if dba is None:
        #             continue


        # If Depth Averaged Velocity (DAV) data available, (i.e. any of the
        # `*_water_vx/vy` sensors are in the data) get the values and calculate
        # the mean position and time of the segment as the postion and time
        # for the velocity estimation
        if file_check.dav_sensors:
            # get segment mean time, segment mean lat, and segment mean lon
            # (underwater portion only)
            seg_time, seg_lat, seg_lon = ooidac.processing.velocity.get_segment_time_and_pos(
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
            vx, vy = ooidac.processing.velocity.get_u_and_v(dba, check_files=next2files)

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
                status.add_nc(out_nc_file, dba_file)

        status.add_src(dba_file)

    # Delete the temporary directory once files have been moved
    try:
        logging.debug('Removing temporary files:')
        shutil.rmtree(ncw.tmp_dir)
    except OSError as e:
        logging.error(e)
        logging.error('manually try to remove temporary directory {:s}'.format(
            ncw.tmp_dir)
        )

    # write the processed files and last profile id to status.json
    logging.debug('Writing run status to status.json')
    if start_profile_id > 0:
        status.info['next_profile_id'] = ncw.profile_id
        status.write()
    # already_processed = set(status['files_processed'])
    # set_processed_dbas = set(processed_dbas)
    # processed_dbas = list(set_processed_dbas.difference(already_processed))
    # processed_dbas.sort(key=sort_function)  # sets don't preserve order
    # status['files_processed'].extend(processed_dbas)
    # status['profiles_created'].extend(output_nc_files)
    # status['profile_to_data_map'].extend(profile_to_data_map)
    # with open(status_path, 'w') as fid:
    #     json.dump(status, fid, indent=2)

    # Print the list of files created
    sys.stdout.write('Profiles NC files written:\n')
    for output_nc_file, source_dba in status.data_map:
        if output_nc_file not in already_processed:
            sys.stdout.write('\t{:s} -> {:s}\n'.format(
                source_dba, output_nc_file))

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
        description=str(main.__doc__),
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
