#!/usr/bin/env python

import os
import sys
import logging
import argparse
import tempfile
import datetime
import shutil
# import pdb
import glob
# temporary addition to test this script -SP 2019-01-30 to be able to load
# gncutils -SP 2019-03-11 Maybe don't need after all.  We'll see.
sys.path.append('C:\\Users\\spearce\\code\\python\\gliderdac\\')
from ooidac.writers.slocum.ooi_ProfileNetCDFWriter import \
    ProfileNetCDFWriter

from ooidac.constants import NETCDF_FORMATS, LLAT_SENSORS
from ooidac.validate import validate_sensors, validate_ngdac_var_names
import numpy as np

import ooidac.processing as process_dba
from ooidac.data_classes import DbaData
from ooidac.profiles import Profiles
from dba_file_sorter import sort_function


def main(args):
    """Parse one or more Slocum glider ascii dba files and write IOOS
    NGDAC-compliant Profile NetCDF files
    """

    # Set up logger
    log_level = getattr(logging, args.loglevel.upper())
    log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
    logging.basicConfig(format=log_format, level=log_level)

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
        # sort dba files since the dbd format doesn't sort well.
        dba_files.sort(key=sort_function)
    else:
        dba_files = args.dba_files

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
    ncw = ProfileNetCDFWriter(config_path, comp_level=comp_level,
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
    tmp_dir = tempfile.mkdtemp()
    logging.debug('Temporary NetCDF directory: {:s}'.format(tmp_dir))

    # Write one NetCDF file for each input file
    output_nc_files = []
    processed_dbas = []
    for dba_file in dba_files:

        if not os.path.isfile(dba_file):
            logging.error('Invalid dba file specified: {:s}'.format(dba_file))
            continue

        logging.info('Processing dba file: {:s}'.format(dba_file))

        # Parse the dba file
        dba = DbaData(dba_file)

        if dba is None or dba.N == 0:
            logging.warning('Skipping empty dba file: {:s}'.format(dba_file))
            continue

        # Check on the "goodness" of the dba data file.  I.e. checks that
        # there is any science data and discovers which timestamp to use and
        # file_good = dba.check_goodness()
        # if not file_good:
        #     logging.warning(
        #         ('The file {:s} does not have the needed '
        #          'variables or data to proceed').format(dba_file)
        #     )
        #     continue

        # %------ Process DBA ------%
        # build the llat (longitude, latitude, altitude, time;
        # i.e. x,y,z,t) variables for NetCDF coordinates
        dba = process_dba.create_llat_sensors(dba)

        # use conductivity, temperature, etc. to create salinity and density
        # variables in the dba object
        dba = process_dba.ctd_data(dba, ctd_sensors, ncw)

        # get the u and v water velocity components (it must come from the
        # next seqment file as it is not processed until as segment is complete)
        # TODO: put code here to get depth averaged water velocity. Also make
        #  code to grab the median lat, lon, and time for the whole segment
        #  to use for u and v. Maybe that code will be in process_dba or the
        #  dba parsing class.

        # %------ Build Profiles ------%
        # Use the `Profiles` class to find the inflection points and split
        # the data into profiles
        profiles = Profiles(dba)

        profiles.find_profiles_by_depth(10)

        # applies the filter functions in profile_filters.py to filter profiles
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

        # TODO: is this done elsewhere? Ans: yes in create_llat_sensors
        # Calculate depth from pressure and replace the old llat_depth
        # zi = [s['sensor_name'] for s in dba['sensors']].index('llat_depth')
        # dba['data'][:, zi] = calculate_depth(ctd_data[:, 1], mean_lat)

        # %------ Write Profiles to NetCDF ------%
        for profile in profiles:
            # TODO: re-write this section in terms of profile indices stored
            profile_times = profile.getdata('llat_time')
            # Calculate and convert profile mean time to a datetime
            prof_start_time = float(profile_times[0])
            mean_profile_epoch = float(np.nanmean([profile_times[0],
                                                   profile_times[-1]]))
            if np.isnan(mean_profile_epoch):
                logging.warning('Profile mean timestamp is Nan')
                continue
            # If no start profile id was specified on the command line,
            # use the mean_profile_epoch as the profile_id since it will be
            # unique to this profile and deployment
            if args.start_profile_id < 1:
                ncw.profile_id = int(mean_profile_epoch)
            pro_mean_dt = datetime.datetime.utcfromtimestamp(mean_profile_epoch)
            prof_start_dt = datetime.datetime.utcfromtimestamp(prof_start_time)

            # TODO: I feel like all of this should be a method of the writer.

            # Create the output NetCDF path
            pro_mean_ts = pro_mean_dt.strftime('%Y%m%dT%H%M%SZ')
            prof_start_ts = prof_start_dt.strftime('%Y%m%dT%H%M%SZ')
            if dba.file_metadata['filename_extension'] == 'dbd':
                filetype = 'delayed'
            elif dba.file_metadata['filename_extension'] == 'sbd':
                filetype = 'rt'
            else:
                logging.warning(
                    'Unknown filename extension {:s}, {:s}'.format(
                        dba.source_file, dba.file_metadata['filename_extension']
                    )
                )
                continue
            profile_filename = '{:s}-{:s}-{:s}'.format(
                ncw.attributes['deployment']['glider'], prof_start_ts,
                filetype)
            # Path to temporarily hold file while we create it
            tmp_fid, tmp_nc = tempfile.mkstemp(
                dir=tmp_dir, suffix='.nc',
                prefix=os.path.basename(profile_filename)
            )
            os.close(tmp_fid)

            out_nc_file = os.path.join(output_path, '{:s}.nc'.format(
                profile_filename))
            if os.path.isfile(out_nc_file):
                if args.clobber:
                    logging.info(
                        'Clobbering existing NetCDF: {:s}'.format(out_nc_file))
                else:
                    logging.warning(
                        'Skipping existing NetCDF: {:s}'.format(out_nc_file))
                    continue

            # Initialize the temporary NetCDF file
            try:
                ncw.init_nc(tmp_nc)
            except (OSError, IOError) as e:
                logging.error('Error initializing {:s}: {}'.format(tmp_nc, e))
                continue

            try:
                ncw.open_nc()
                # Add command line call used to create the file
                ncw.update_history('{:s} {:s}'.format(sys.argv[0], dba_file))
            except (OSError, IOError) as e:
                logging.error('Error opening {:s}: {}'.format(tmp_nc, e))
                os.unlink(tmp_nc)
                continue

            # Create and set the trajectory
            # trajectory_string = '{:s}'.format(ncw.trajectory)
            ncw.set_trajectory_id()
            # Update the global title attribute with the name of the source
            # dba file
            ncw.set_title(
                '{:s}-{:s} Vertical Profile'.format(
                    ncw.deployment_configs['glider'],
                    pro_mean_ts
                )
            )

            # Create the source file scalar variable
            ncw.set_source_file_var(
                dba.file_metadata['filename_label'], dba.file_metadata)

            # Update the self.nc_sensors_defs with the dba sensor definitions
            ncw.update_data_file_sensor_defs(dba['sensors'])

            # Find and set container variables
            ncw.set_container_variables()

            # Create variables and add data
            for var_name in dba.sensor_names:
                var_data = dba[var_name]['data'][profile_ii]
                logging.debug('Inserting {:s} data array'.format(var_name))
                ncw.insert_var_data(var_name, var_data)

            # Write scalar profile variable and permanently close the NetCDF
            # file
            nc_file = ncw.finish_nc()

            if nc_file:
                try:
                    shutil.move(tmp_nc, out_nc_file)
                    os.chmod(out_nc_file, 0o755)
                except IOError as e:
                    logging.error(
                        'Error moving temp NetCDF file {:s}: {:}'.format(
                            tmp_nc, e)
                    )
                    continue

            output_nc_files.append(out_nc_file)

        processed_dbas.append(dba_file)

    # Delete the temporary directory once files have been moved
    try:
        logging.debug('Removing temporary directory: {:s}'.format(tmp_dir))
        shutil.rmtree(tmp_dir)
    except OSError as e:
        logging.error(e)
        return 1

    # Print the list of files created
    for output_nc_file in output_nc_files:
        os.chmod(output_nc_file, 0o664)
        sys.stdout.write('{:s}\n'.format(output_nc_file))

    return 0


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
