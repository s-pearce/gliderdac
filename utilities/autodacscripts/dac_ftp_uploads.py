# Glider DAC FTP upload

import os
import sys
import logging
import json
import glob
import netrc
from functools import partial
import traceback
import ftplib
from ftplib import FTP
import numpy as np
from report_email import ReportEmail

"""dac_ftp_uploads.py
Upload profile netCDF files to the GliderDAC via FTP

Used with `gliderdac` code to upload processed profile netCDF files over an FTP
connection to the IOOS National GliderDAC.  Used in conjunction with the 
gliderdac created status.json file in the associated configuration directory, 
`dac_ftp_uploads` will transfer processed profiles that have not yet been
uploaded to the DAC.  The ability to upload only un-uploaded profiles allows
for near real-time, newly produced data from currently deployed gliders to be
pushed to the DAC regularly.

"""
__version__ = "$Revision: 0.2.0 $"
__author__ = "Stuart Pearce (stuart.pearce@oregonstate.edu)"
__date__ = "$Date: 2020-03-05 15:09 $"


# ToDo: import configuration information from config file
# ToDo: make option arguments so that the debug level can be selected
# ToDo: check nc_list with status.json profiles_uploaded to avoid unneccessary
#       FTP connections
# ToDo: add logger to the report_email module too.
# ToDo: Cross-check using status.json to compare nc files to "profiles_created".

USERNAME, _, DACPASSWD = netrc.netrc().authenticators('data.ioos.us')

# Set up logging
log_level = getattr(logging, 'DEBUG')
log_format = '%(module)s:%(levelname)s:%(message)s'
logging.basicConfig(format=log_format, level=log_level)

# emailer for sending out error reports
ereport = ReportEmail()

def traceback_printer():
    """A shorter error report than the typical error raise
    
    returns
    -------
    err_str : str
        a string with the error, error text and the traceback information for
        an error.
    """
    etype, val, tb = sys.exc_info()
    xtb = traceback.extract_tb(tb)
    fexc = traceback.format_exception(etype, val, tb)
    err_str = "{:s}".format(fexc[-1])
    err_str += "line:{:d} in {:s}: {:s}\n".format(
        xtb[0].lineno, xtb[0].name, xtb[0].line)
    if len(xtb) > 1:
        err_str += "  traceback:\n"
        for ii in range(1,len(xtb)):
            err_str += "    line {:d} {:s}: {:s}".format(
                xtb[ii].lineno, xtb[ii].name, xtb[ii].line)
    return err_str

# 1. Receive Inputs config directory and directory of deployment profile netCDFs

# inputs will be configuration directory and netCDF directory
config_dir = sys.argv[1]
nc_dir = sys.argv[2]  # or do we want a list or glob?
logging.debug("config_dir: {:s}\nnc_dir: {:s}".format(config_dir, nc_dir))

# check that the directories exist or else report out they are missing
if not os.path.exists(config_dir) or not os.path.exists(nc_dir):
    # exit and report out no directories found
    ereport.sendmail(
        "Either the config or nc directory does not exist for an attempted "
        "ftp GliderDAC upload.\n"
        "Config directory: {:s}\nnc directory: {:s}".format(config_dir, nc_dir)
    )
    sys.exit(1)

nc_dir = os.path.realpath(nc_dir)

# 2. check if new profile nc files exist for glider-deployment
status_json_file = os.path.join(config_dir, 'status.json')
logging.debug("status.json: {:s}".format(status_json_file))

if not os.path.exists(status_json_file):
    # report out a missing status.json and exit
    ereport.sendmail(
        "The status.json file is missing for an attempted ftp GliderDAC upload"
        "\n {:s}".format(status_json_file)
    )
    sys.exit(1)

with open(status_json_file, 'r') as fid:
    status = json.load(fid)
    logging.debug("status.json opened")

profiles_uploaded = status['profiles_uploaded']
logging.debug("len(profiles_uploaded): {:d}".format(len(profiles_uploaded)))

nc_glob = os.path.join(nc_dir, '*.nc')
local_nc_list = glob.glob(nc_glob)
logging.debug("len(local_nc_list): {:d}".format(len(local_nc_list)))

if not local_nc_list:
    # report out there are no nc files found and exit
    ereport.sendmail(
        "There does not appear to be any NC files available to upload from the"
        "directory:\n {:s}".format(nc_dir)
    )
    sys.exit(1)

# Just want the list to be the basenames for comparison to files on the DAC
local_nc_list = list(map(os.path.basename, local_nc_list))

# 3. if new nc files exist get the DAC directory name (trajectory name)
# Note: this could just call the existing script return_deployment.py but
#       they would both need to be implemented into the GliderDAC code.

if local_nc_list == profiles_uploaded:
    logging.info("No new profiles to upload. Exiting.")
    sys.exit(0)

deployment_file = os.path.join(config_dir, 'deployment.json')
with open(deployment_file, 'r') as fid:
    deploy_info = json.load(fid)
dac_traj = deploy_info['trajectory_name']  
# need to add "-delayed" to the trajectory name
try:
    mode = deploy_info['global_attributes']['mode']
except KeyError as e:
    tracestr = traceback_printer()
    error_msg = (
        "The deployment.json file {:s} does not have a 'mode' attribute.\n"
        "{:s}".format(deployment_file, tracestr)
    )
    logging.warning(error_msg)
    ereport.sendmail(error_msg)
    raise(e)
    
if mode == 'rt':
    dac_dir = dac_traj
elif mode == 'delayed':
    dac_dir = dac_traj + "-" + mode
else:
    # report mode not available and quit
    ereport.sendmail(
        "The deployment.json file does not have a correct 'mode' attribute:\n"
        "{:s}".format(deployment_file)
    )
    sys.exit(1)
logging.debug("trajectory directory: {:s}".format(dac_dir))

# 4. start FTP connection
logging.debug("attempting FTP connection...")
try:
    ftp = FTP('data.ioos.us', user=USERNAME, passwd=DACPASSWD)
except ftplib.all_errors as e:
    # report out connection cannot establish and raise error
    tracestr = traceback_printer()
    error_msg = (
        "FTP connection failed for attempted GliderDAC upload for nc files"
        "from {:s}\n{:s}".format(nc_dir, tracestr)
    )
    logging.error(error_msg)
    ereport.sendmail(error_msg)
    raise(e)
finally:
    del DACPASSWD

with ftp:
    # 5. Check that the directory exists.
    # Note: use a try except statement and catch error `error_perm` if
    #       the directory does not exist
    # either:
    try:
        resp = ftp.cwd(dac_dir)
    except ftplib.all_errors as e:
        # report out DAC trajectory not created/ directory DNE
        dirlist = ftp.nlst()
        tracestr = traceback_printer()
        if dac_dir in dirlist:
            # report out DAC ftp directory change failure
            error_msg = (
                "During an attempted GliderDAC FTP upload, the directory {:s}"
                "exists, but FTP was unable to change to that "
                "directory.\n{:s}".format(dac_dir, tracestr)
            )
            logging.error(error_msg)
            ereport.sendmail(error_msg)
            raise(e)
        else:
            # report out DAC trajectory not created/ directory DNE
            error_msg = (
                "The GliderDAC trajectory {:s} has not been created or does "
                "not exist\n{:s}".format(dac_dir, tracestr)
            )
            logging.error(error_msg)
            ereport.sendmail(error_msg)
            raise(e)

    logging.debug("successfully changed directory")

    # S?. Retrieve list of files already uploaded and take the difference
    #     between remote files and local files to get the files needing upload
    try:
        dac_nc_list = ftp.nlst('*.nc')
    except ftplib.all_errors as e:
        # report out the error and exit
        tracestr = traceback_printer()
        error_msg = (
            "The GliderDAC trajectory {:s} has failed at giving a netCDF file "
            "listing\n {:s}".format(dac_dir, tracestr)
        )
        logging.warning(error_msg)
        ereport.sendmail(error_msg)
        raise(e)

    logging.debug("received dac_nc_list, len: {:d}".format(len(dac_nc_list)))
    
    # Find the nc files that have not been uploaded to the DAC by comparison
    # with the local nc file list.
    
    # The partial function here partially fills out os.path.join with the nc_dir
    # kind of like a conditional default to shorten the input for use in the map
    # function below
    #ncjoin = partial(os.path.join, nc_dir)
    # now even though these dac_nc_list is from the remote directory on the DAC, 
    # this creates a fake full path list so that I can compare the local_nc_list
    # with the dac_nc_list and the remainder from the set difference will be
    # the files that need to be uploaded already in full path format.  This is
    # better than comparing the basefile names and then reconstructing the full
    # path name from the difference.
    #fake_full_path_dac_list = list(map(ncjoin, dac_nc_list))
    files_to_upload = np.setdiff1d(local_nc_list, dac_nc_list)
    logging.debug(
        "first estimation of number of files to upload: {:d}".format(
            len(files_to_upload)))
    
    # Check that dac_nc_list is the same as the status profiles_uploaded
    # list, if not, add the missing files to the upload list

    # For error handling, we want to compare the kept profiles_uploaded list
    # with the remote directory listing from the DAC.  If there is a file in 
    # one that is not in the other, that file should be uploaded/re-uploaded.
    # This should take care of the 2 cases where the ftp connection may have
    # been severed during an upload; where some portion of a file has been
    # uploaded to the DAC, but is incomplete/corrupted, or where this program
    # thinks it has uploaded correctly, but the file is not found in the remote
    # DAC directory.
    files_to_reupload = np.setxor1d(profiles_uploaded , dac_nc_list)
    if len(files_to_reupload) > 0:
        for file in files_to_reupload:
            logging.debug("{:s} needs to be re-uploaded".format(file))
    files_to_upload = np.append(files_to_upload, files_to_reupload)
    logging.debug(
        "second estimate of num of files to upload: {:d}".format(
            len(files_to_upload)))

    # 7. batch upload new profile nc files in binary mode
    for file in files_to_upload:
        #base_filename = os.path.basename(file)
        filepath = os.path.join(nc_dir, file)
        logging.debug("Attempting to upload file: {:s}".format(file))
        try:
            with open(filepath, 'rb') as fid:
                resp = ftp.storbinary('STOR {:s}'.format(file), fid)
        except ftplib.all_errors as e:
            tracestr = traceback_printer()
            error_msg = (
                "Attempted GliderDAC file {:s} FTP upload failed:\n"
                "{:s}".format(file, tracestr)
            )
            logging.error(error_msg)
            ereport.sendmail(error_msg)
        except OSError as e:
            tracestr = traceback_printer()
            error_msg = (
                "Attempted GliderDAC file {:s} FTP upload failed because file "
                "could not be opened:\n{:s}".format(file, tracestr)
            )
            logging.error(error_msg)
            ereport.sendmail(error_msg)
        else:
            if not resp == '226 Transfer complete.':
                try:
                    ftp.delete(file)
                except ftplib.all_errors as e:
                    tracestr = traceback_printer()
                    error_msg = (
                        "Attempt to delete uploaded GliderDAC file {:s} failed.\n"
                        "{:s}".format(file, tracestr)
                    )
                    logging.error(error_msg)
                    ereport.sendmail(error_msg)
            else:
                logging.info(
                    'Uploaded {:s}'.format(file))
                status['profiles_uploaded'].append(file)
            

# 8. catalog the profiles that have been uploaded so they aren't reuploaded
logging.debug("Writing out new status.json file")
with open(status_json_file, 'w') as fid:
    json.dump(status, fid)
