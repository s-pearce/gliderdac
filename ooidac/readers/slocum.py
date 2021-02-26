"""Routines for parsing Slocum glider ascii dba files created with the
shoreside executables provided by Teledyne Webb Research.  These routines will
not parse Slocum matlab (.m,.dat) files produced by dba2_orig_matlab or
dba2_glider_data"""

import os
import logging
import numpy as np
import time

logger = logging.getLogger(os.path.basename(__name__))


# ToDo: bring comments up to date if necessary
def parse_dba(dba_file, fast=False):
    """Parse a Slocum dba ascii table file.

    Args:
        dba_file: dba file to parse

    Returns: A dictionary containing the file metadata, sensor defintions
    and data
    """
    if not os.path.isfile(dba_file):
        logging.error('Invalid dba file: {:s}'.format(dba_file))
        return
    t0 = time.time()
    try:
        with open(dba_file, 'r') as dbafid:
            # Parse the dba header
            dba_headers = _parse_dba_header(dbafid)
            # Parse the dba sensor definitions
            sensors, sensor_defs = _parse_dba_sensor_defs(dbafid)

    except IOError as e:
        logging.error('Error opening {:s} dba file: {}'.format(
            dba_file, e))
        return

    if not dba_headers or not sensor_defs:
        return

    # Figure out what line the data table begins on using the header's
    # num_ascii_tags + num_lable_lines
    if 'num_label_lines' not in dba_headers:
        logger.warning('num_label_lines header missing: {:s}'.format(
            dba_file))
        return
    num_header_lines = int(dba_headers['num_ascii_tags'])
    num_label_lines = int(dba_headers['num_label_lines'])
    num_columns = int(dba_headers['sensors_per_cycle'])
    # Total number of header lines before the data matrix starts
    total_header_lines = num_header_lines + num_label_lines

    # Parse the ascii table portion of the dba file
    if fast:
        data = _fast_load_dba_data(dba_file, total_header_lines)
    else:
        data = _load_dba_data(dba_file, total_header_lines)

    if data is None or len(data) == 0:
        logger.warning('Data length is 0 in dba file: {:s}'.format(
            dba_file))
    elif num_columns != data.shape[1]:
        logger.warning(
            'Glider data file does not have the same'
            'number of columns as described in header.\n'
            'described {:d}, actual {:d}'.format(
                num_columns, data.shape[1])
        )
    t1 = time.time()
    logger.debug("Time elapsed for parser, {:0.0f}".format(t1 - t0))

    dba = {'header': dba_headers, 'sensor_names': sensors,
           'sensor_defs': sensor_defs, 'data': data}
    return dba


def _parse_dba_header(fid):
    """Parse the header information in a Slocum dba ascii table file.
    All header lines of the format 'key: value' are parsed.

    Args:
        fid: dba file object to parse

    Returns:
        A dictionary containing the file metadata
    """
    dba_headers = {}
    file_byte_pos = fid.tell()
    try:
        line = fid.readline()
        lines_read = 1
        while line:
            tokens = line.strip().split(': ')
            if len(tokens) != 2:
                fid.seek(file_byte_pos)
                lines_read -= 1
                break

            dba_headers[tokens[0]] = tokens[1]
            file_byte_pos = fid.tell()
            line = fid.readline()
            lines_read += 1
    except IOError as e:
        logging.error('Error parsing {:s} dba header: {}'.format(
            fid.name, e))
        return

    if not dba_headers:
        logging.warning('No headers parsed: {:s}'.format(
            fid.name))
        return

    if 'num_ascii_tags' in dba_headers:
        if lines_read != int(dba_headers['num_ascii_tags']):
            logging.warning(
                'Unexpected number of header fields: {:s}'.format(fid.name))
            return
    else:
        logging.warning('num_ascii_tags header line missing: {:s}'.format(
            fid.name))
        return

    # Add the full path to the dba_file
    dba_headers['full_path'] = os.path.realpath(fid.name)
    dba_headers['source_file'] = os.path.basename(fid.name)
    # Add the dba file size
    dba_headers['file_size_bytes'] = os.stat(fid.name).st_size

    return dba_headers


def _parse_dba_sensor_defs(fid):
    """Parse the sensor definitions in a Slocum dba ascii table file.

    Args:
        fid: dba file object to parse

    Returns:
        An array of dictionaries containing the file sensor definitions
    """

    try:
        # Get the sensor names line
        sensors_line = fid.readline().strip()
        # Get the sensor units line
        units_line = fid.readline().strip()
        # Get the datatype byte storage information
        bytes_line = fid.readline().strip()
    except IOError as e:
        logging.error(
            'Error parsing {:s} dba header: {:s}'.format(fid.name, e))
        return None, None

    sensors = sensors_line.split()
    units = units_line.split()
    datatype_bytes = bytes_line.split()

    if not sensors:
        logging.warning(
            'No sensor defintions parsed: {:s}'.format(fid.name))
        return None, None

    sensor_defs = {}
    for ii in range(len(sensors)):
        sensor_defs[sensors[ii]] = {
            'sensor_name': sensors[ii],
            'attrs': {
                'units': units[ii], 'bytes': int(datatype_bytes[ii]),
                'source_sensor': sensors[ii], 'long_name': sensors[ii]
            },
        }
    return sensors, sensor_defs


def _load_dba_data(dba_file, num_header_lines=17):

    # Use numpy.loadtxt to load the ascii table, skipping header rows and
    # requiring a 2-D output array
    try:
        t0 = time.time()
        data_table = np.loadtxt(dba_file, skiprows=num_header_lines,
                                ndmin=2)
        t1 = time.time()
        elapsed_time = t1 - t0
        logger.debug('DBD parsed in {:0.0f} seconds'.format(
            elapsed_time))
    except ValueError as e:
        logger.warning('Error parsing {:s} ascii data table: {:s}'.format(
            dba_file, e))
        return

    return data_table


def _fast_load_dba_data(dba_file, num_header_lines=17):

    # Use numpy.loadtxt to load the ascii table, skipping header rows and
    # requiring a 2-D output array
    try:
        t0 = time.time()
        # read each row of data & use np.array's ability to grab a
        # column of an array
        data = []
        # pdb.set_trace()
        line_num = 1
        with open(dba_file, 'r') as dba_fid:
            for line in dba_fid.readlines():
                if line_num > num_header_lines:
                    data.append(line.split())
                line_num += 1
        t1 = time.time()
        elapsed_time = t1 - t0
        logger.debug('DBD parsed in {:0.0f} seconds'.format(
            elapsed_time))
        data_array = np.array(
            data, dtype=np.float)  # NOTE: this is an array of strings
    except ValueError as e:
        logger.warning('Error parsing {:s} ascii data table: {:s}'.format(
            dba_file, e))
        return

    return data_array


def parse_dba_header(dba_file):
    if not os.path.isfile(dba_file):
        logging.error('Invalid dba file: {:s}'.format(dba_file))
        return
    with open(dba_file, 'r') as dba_fid:
        dba_header = _parse_dba_header(dba_fid)
    return dba_header


def parse_dba_sensor_defs(dba_file):
    if not os.path.isfile(dba_file):
        logging.error('Invalid dba file: {:s}'.format(dba_file))
        return
    with open(dba_file, 'r') as dba_fid:
        dba_header = _parse_dba_header(dba_fid)
        dba_header['sensors'] = _parse_dba_sensor_defs(dba_fid)
    return dba_header
