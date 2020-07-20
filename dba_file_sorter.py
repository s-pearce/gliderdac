import os
import re
import logging

logger = logging.getLogger(os.path.basename(__name__))

slocumregex = (
    r'[a-z_0-9]+(?:_|-)(\d{4}(?:_|-)\d{3}(?:_|-))(\d{1,2})(?:_|-)'
    r'(\d{1,4})\..+')
regex = re.compile(slocumregex)


def sort_function(filename):
    fname = os.path.basename(filename)
    match = regex.search(fname)
    if not match:
        error_msg = (
            'File {:s} does not match the Slocum filename format.'.format(
                filename
            ))
        logging.warning(error_msg)
        raise ValueError(error_msg)
    else:
        mission_num = '{:02d}'.format(int(match.group(2)))
        segment_num = '{:04d}'.format(int(match.group(3)))
        sort_str = match.group(1) + mission_num + '_' + segment_num
    return sort_str
