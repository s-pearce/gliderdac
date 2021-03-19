#! 
"""DAC NC utilities dac_nc_utils.py

a module of utilities for post-processed NGDAC netCDF files for update or
correction.

Author: Stuart Pearce
Email: stuart.pearce at oregonstate.edu
"""

import os
import glob
import re
import datetime as dt
import numpy as np
import netCDF4

__author__ = "Stuart Pearce"
__email__ = "stuart.pearce@oregonstate.edu"
__date__ = "2020-05-12"
__version__ = 0.4

# ToDo: make a separate function that handles flag values based on text or
#       integer input and checks that all are valid

# The FLAG_DICT is created this way only because the full FLAG_MEANINGS
# and FLAG_VALUES are inputs to a QC variable attributes and reduces
# repeated typing.
FLAG_MEANINGS = [
            "no_qc_performed", 
            "good_data", 
            "probably_good_data",
            "bad_data_that_are_potentially_correctable",
            "bad_data",
            "value_changed",
            "not_used",
            "not_used",
            "interpolated_value",
            "missing_value"]
FLAG_VALUES = [0,1,2,3,4,5,6,7,8,9]

FLAG_DICT = dict(zip(FLAG_MEANINGS, FLAG_VALUES))
# ^ drops one of the 'not_used' because of redundant keywords.
FLAG_DICT.pop('not_used')  # removes 'not_used' entirely so it is not returned


def handle_file_type(func):
    """decorator function for netCDF file object.
    Takes the fileobj input for any of these netCDF file operations, and if the
    input is a string filename, opens the netCDF file in append mode.  If 
    The input is a netCDF object already, it just passes it to the utility
    function.
    """
    def wrapper(fileobj, *args, **kwargs):
        if isinstance(fileobj, str):
            if not os.path.exists(fileobj):
                raise FileNotFoundError(
                    "File {:s} does not exist".format(file))
            with netCDF4.Dataset(fileobj, 'a') as nc:
                func(nc, *args, **kwargs)
        if isinstance(fileobj, netCDF4._netCDF4.Dataset):
            #_assert_nc_write_mode(fileobj)  # excluding for now
            func(fileobj, *args, **kwargs)
    wrapper.__doc__ = func.__doc__
    return wrapper


@handle_file_type
def flag_data(ncfile, varname, flag_value, start=None, stop=None):
    """
    INCOMPLETE
    
    params
    ------
    ncfile : string or netCDF4.Dataset instance open for appending
        The netCDF file to be appended
    varname : string
        The variable name that is to be flagged
    flag_value : string or int
        A valid NGDAC qc variable flag value
    start : datetime object
        The starting date and time of the affected data that the flag value
        describes
    stop : datetime object
        The ending date and time of the affected data that the flag value
        describes
    

    returns
    -------
    None : makes change to ncfile, if ncfile is a string filename, the file is
           opened, appended, and closed.
    """
    # ToDo: make a more overall function something like `flag_data` with inputs
    #       file, variable, flag_value, and optional daterange.  Based on flag
    #       value it will fill the variable or not and create the qc flag variable
    #       for the time range if specified. adds a comment and updates history
    #       all together.

    raise NotImplementedError("this function has not been completed yet.")

    if isinstance(flag_value, str) and flag_value in FLAG_DICT:
        flag = FLAG_DICT[flag_value]
    elif isinstance(flag_value, int) and flag_value in FLAG_DICT.values():
        flag = flag_value
    else:
        raise ValueError(
            "flag_value must be a valid flag_value keyword or integer value")

    if start is None:
        pass
    if stop is None:
        pass
    
    if flag == 4:
        fill_variable(ncfile, varname)
    add_qc_var(ncfile, varname, flag)
    update_history(ncfile, history_comment)
    update_comment(ncfile, varname, comment)


@handle_file_type
def fill_variable(ncfile, varname, fill_value=None):
    """Fill a variable with a fill value
    
    `fill_variable` will rewrite the data in a netCDF variable, filling with
    the fill value (either provided by `fill_value`, using already established
    fill values or using the default fill value for the data type of the
    variable.
    https://ioos.github.io/ioosngdac/ngdac-netcdf-file-format-version-2.html
    
    If 'fill_value' is given as a single value, `add_qc_var` will fill
    out the QC flag variable to the same length as the `varname` data.
    
    Note: Since this module applies only to Glider DAC profile netCDFs,
    seeing how it will rarely require a fill to apply to less than a full
    profile, at this time, `fill_variable` only fills the entire profile netCDF
    file.  It may require change at a later time.
    
    params
    ------
    ncfile : string or netCDF.Dataset instance open in append mode
        netCDF file to append
    varname : string.
        The name of the variable to fill with fill values
    fill_value : (optional) int, float, or nan. 
        The fill value number used to fill the variable

    
    returns
    -------
    None : If ncfile is a string filename, the file is written and closed
    """
    # ToDo: add datetime range
    
    # file opening error checking is handled by the handle_file_type decorator
    # and provides this function with an open for appending netCDF file object
    
    # Error check for the netCDF variable
    if varname not in ncfile.variables:
        raise IndexError(
            "Variable {:s} not found in NetCDF file"
            "{:s}".format(varname, ncfile.file_path()))
    
    # If a new fill value is not given, check the fill value already used. 
    # If none was used, use the typical fill value for the data type.
    if fill_value:
        fill_val = fill_value
    elif "_FillValue" in ncfile[varname].ncattrs():  # newer attribute name
        fill_val = ncfile[varname]._FillValue
    elif "missing_value" in ncfile[varname].ncattrs():  # older attribute name
        fill_val = ncfile[varname].missing_value
    else:
        fill_val = netCDF4.default_fillvals[ncfile[varname].dtype.str[1:]]
    
    # replace the data with fill values
    size = ncfile[varname].size
    
    # no `try` statement so the natural RuntimeError will be raised if the file
    # is not open in append mode
    ncfile[varname][:] = np.full(size, fill_val)


@handle_file_type
def add_qc_var(ncfile, varname, qc_flags):
    """Add a QC variable using provided flags
    
    `add_qc_var` creates and adds a QC flag variable that describes the 
    variable `varname` to Glider DAC netCDF file `ncfile` using provided 
    `qc_flags`.  The format of the QC variable follows the format examples 
    described in IOOS_Glider_netCDF version 2.0.
    https://ioos.github.io/ioosngdac/ngdac-netcdf-file-format-version-2.html
    
    If `qc_flags` is given as a single value, `add_qc_var` will fill
    out the QC flag variable to the same length as the `varname` data.
    
    DAC QC Flags table of values:
    #       Flag
    -       ----
    0	    no_qc_performed
    1	    good_data
    2	    probably_good_data
    3	    bad_data_that_are_potentially_correctable
    4	    bad_data
    5	    value_changed
    6	    not_used
    7	    not_used
    8	    interpolated_value
    9	    missing_value
    -127    netCDF fill value

    
    params
    ------
    ncfile : netCDF file to append
    varname : string, the variable name that the qc flag variable represents.
        The QC variable name will be `varname` suffixed with "_qc" at the end.
    qc_flags : array-like, an array the same length as the variable `varname`
        or a scalar value to fill the qc variable.
    
    returns
    -------
    None : netCDF file is written back to the file.
    """
    # ToDo: add a datetime range
    # ToDo: check that varname_qc does not already exist and check the flags on
    #       it.  If nothing to do then return
    
    # file opening error checking is handled by the handle_file_type decorator
    # and provides this function with an open for appending netCDF file object
    
    # Check that qc_flags fits the DAC QC flag description
    qc_flags = np.atleast_1d(qc_flags)
    if (
            np.any(qc_flags>9) or np.any(qc_flags<-127) or
            np.any(np.logical_and(qc_flags<0, qc_flags>-127))):
        raise ValueError(
            "qc_flags can only be values 0-9 or a fill value of -127")
    
    # A couple of error checks for the netCDF variable
    
    if varname not in ncfile.variables:
        raise IndexError(
            "Variable {:s} not found in NetCDF file"
            "{:s}".format(varname, ncfile.file_path()))
    if "standard_name" not in ncfile[varname].ncattrs():
        raise KeyError("'standard_name' is a required variable attribute")
    if "long_name" not in ncfile[varname].ncattrs():
        raise KeyError("'long_name' is a required variable attribute")
    if "ancillary_variables" not in ncfile[varname].ncattrs():
        print("'ancillary_variables' is a required variable attribute. Adding")
        # no `try` statement here so the natural AttributeError will be raised
        # if the file is not open in append mode
        ncfile[varname].setncattr('ancillary_variables', '')

    # get attributes and characteristics for use in the new QC variable
    dims = ncfile[varname].dimensions
    long_name = ncfile[varname].long_name
    stand_name = ncfile[varname].standard_name
    filters = ncfile[varname].filters()
    zlib_ = filters['zlib']
    comp = filters['complevel']
    size = ncfile[varname].size

    # check qc_flags input
    if len(qc_flags) == 1:
        qc_flags = np.full(size, qc_flags)
    if len(qc_flags) != size:
        raise ValueError(
            "qc_flags is not a scalar or the same length as variable "
            "{:s}".format(varname))

    # create qc variable:
    qc_varname = varname + "_qc"
    if qc_varname not in ncfile.variables:
        
        # no `try` statement so the natural RuntimeError will be raised if the
        # file is not open in append mode
        qc_var = ncfile.createVariable(
            qc_varname, np.byte, dimensions=dims, zlib=zlib_,
            complevel=comp, fill_value=-127)
        
        # write data to it
        qc_var[:] = qc_flags
        
        # write attributes according to DAC format
        qc_var.setncatts({
            "flag_meanings": " ".join(FLAG_MEANINGS),
            "flag_values": FLAG_VALUES,
            "long_name": long_name + " Quality Flag",
            "standard_name": stand_name + " status_flag",
            "valid_max": 9,
            "valid_min": 0
            })
        
        # add the newly created QC variable name to the ancillary_variables
        # attribute in the main variable
        if qc_varname not in ncfile[varname].ancillary_variables:
            ncfile[varname].ancillary_variables += " " + qc_varname

@handle_file_type
def change_global_attr(ncfile, attr, new_value):
    """
    Warning, will overwrite existing attributes, make sure you want to before 
    running.
    """
    # file opening error checking is handled by the handle_file_type decorator
    # and provides this function with an open for appending netCDF file object
    
    # no `try` statement so the natural AttributeError will be raised if the
    # file is not open in append mode
    ncfile.setncattr(attr, new_value)


@handle_file_type
def change_var_attr(ncfile, varname, attr, new_value):
    """
    Warning, will overwrite existing attributes, make sure you want to before
    running.
    """
    # file opening error checking is handled by the handle_file_type decorator
    # and provides this function with an open for appending netCDF file object
    
    # Error check for the netCDF variable
    if varname not in ncfile.variables:
        raise IndexError(
            "Variable {:s} not found in NetCDF file"
            "{:s}".format(varname, ncfile))
    var = ncfile[varname]
    # no `try` statement so the natural AttributeError will be raised if the
    # file is not open in append mode
    var.setncattr(attr, new_value)

    ## uncomment the section below if you ever want the old value to be
    ## a certain value first before commiting to changing the value.  Also
    ## add to the inputs an optional `check_old_value` parameter and remove
    ## the last line above
    # # check if the attr exists first
    # if attr in var.ncattrs():
        # attr_exists = True
    # else:
        # attr_exists = False
    # if (
            # (attr_exists and check_old_value==var.getncattr(attr)) or
            # not attr_exists or
            # check_old_value is None):
        # var.setncattr(attr, new_value)
                


@handle_file_type
def append_var_attr(ncfile, varname, attr, new_value, delim="\n"):
    """Appends a variable's string attribute with additional text
    
    Will concatenate the `new_value` with the attribute's previous value
    delimited by the `delim` character [default is the line feed character \n].
    
    params
    ------
    ncfile : string or netCDF4.Dataset instance
        The netCDF file to update
    varname : string
        The variable name in `ncfile` to update
    attr : string
        The attribute in `varname` to append/update
    new_value : string
        The message/value to append to the current attribute's message/value
    delim : optional string
        The delimiter value to separate the current attribute value with the 
        new supplied attribute value.  Default is the line feed character "\n"
    """
    # file opening error checking is handled by the handle_file_type decorator
    # and provides this function with an open for appending netCDF file object
    
    # Error check for the netCDF variable
    if varname not in ncfile.variables:
        raise IndexError(
            "Variable {:s} not found in NetCDF file"
            "{:s}".format(varname, ncfile.file_path()))
    var = ncfile[varname]
    if attr in var.ncattrs():
        old_value = var.getncattr(attr)
    else:
        old_value = ""
        delim = ""
    value = old_value + delim + new_value
    # no `try` statement here so the natural AttributeError will be raised if
    # the file is not open in append mode
    var.setncattr(attr, value)


@handle_file_type
def update_history(ncfile, message, timestamp=None):
    """
    """
    # file opening error checking is handled by the handle_file_type decorator
    # and provides this function with an open for appending netCDF file object
    
    hist_ts_fmt = '%Y-%m-%dT%H:%M:%SZ'
    
    if timestamp is None:
        now = dt.datetime.utcnow()
        hist_ts = now.strftime(hist_ts_fmt)
    elif isinstance(timestamp, dt.datetime):
        hist_ts = timestamp.strftime(hist_ts_fmt)
    else:
        raise ValueError("timestamp should be a datetime.datetime instance")
    
    if 'history' not in ncfile.ncattrs():
        # no `try` statement so the natural AttributeError will be raised if
        # the file is not open in append mode
        ncfile.setncattr('history', '')
    
    new_history = ncfile.history + '\n' + hist_ts + ' ' + message
    # no `try` statement so the natural AttributeError will be raised if the
    # file is not open in append mode
    ncfile.history = new_history


@handle_file_type
def update_comment(ncfile, message, varname=None):
    """Note: works for GA comments and variable comments?
    """
    # file opening error checking is handled by the handle_file_type decorator
    # and provides this function with an open for appending netCDF file object
    
    if varname is None:
        nclevel = ncfile
    elif varname in ncfile.variables:
        nclevel = ncfile[varname]
    else:
        raise IndexError(
            "Variable {:s} not found in NetCDF file"
            "{:s}".format(varname, ncfile.file_path()))

    if 'comment' not in nclevel.ncattrs():
        nclevel.setncattr('comment', '')

    if nclevel.comment == '':
        prefix = ''
    else:
        prefix = '\n'

    nclevel.comment += prefix + message


regex = re.compile(r'.+_(\d{8}T\d{6}Z)_.+')

def _convert_filename2time(file):
    file = os.path.basename(file)
    match = regex.match(file)
    if match:
        timestr = match.group(1)
        file_dt = dt.datetime.strptime(timestr, '%Y%m%dT%H%M%SZ')
        return file_dt

def select_files_by_daterange(files, start, end=None):
    """
    """
    # ToDo: add the error checks for datetime objects.  Do I want the ability
    #       to have start or end be a string?  I would either need to enforce a
    #       format, or use pandas.

    #if not isinstance(start, dt.datetime):
        
    #if not isinstance(end, dt.datetime):
    
    files = np.array(files)
    file_times = np.array(list(map(_convert_filename2time, files)))
    if end is None:
        end = file_times[-1]
    daterange_indices = np.flatnonzero(np.logical_and(
        file_times >= start,
        file_times <= end))

    # add the file just before the first included file as it might contain the
    # `start` timestamp or timestamps that come after `start`
    first = daterange_indices[0]
    if first > 0:  # only do this if first is not index 0
        first = first-1
        nc = netCDF4.Dataset(files[first])
        with nc:
            ncend = dt.datetime.strptime(
                nc.time_coverage_end, '%Y-%m-%dT%H:%M:%SZ')
            if ncend > start:
                daterange_indices = np.append(
                    np.array([first]), daterange_indices)
    
    return list(files[daterange_indices])


def _assert_nc_write_mode(ncobj):
    if not isinstance(ncobj, netCDF4._netCDF4.Dataset):
        raise TypeError('Input must be an open netCDF4 Dataset instance')
    try:
        ncobj.setncattr('xz_write_test', 'test')
    except AttributeError:
        raise NetCDFOpenModeError(
            "NetCDF object must be open in append mode; 'a' or 'r+'")


class NetCDFOpenModeError(Exception):
    pass


@handle_file_type
def test_dec(ncfile, varname, attr):
    """just a function that tests the handle_file_type decorator function
    """
    var = ncfile[varname]
    if attr not in var.ncattrs():
        print("attribute {:s} not in variable {:s}".format(attr, varname))
    else:
        value = var.getncattr(attr)
        print(
            "Attribute {:s} in variable {:s} value is:\n{}".format(
                attr, varname, value))
