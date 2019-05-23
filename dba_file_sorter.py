import os
import re
import logging

import glob

logger = logging.getLogger(os.path.basename(__name__))


# regex = r'(
#   [_A-Za-z0-9]+)(_|-)(\d{4})(_|-)(\d{3})(_|-)(\d{1,2})(_\-)(\d{1,3})\.'
ooiregex = r'[a-z]{2}_\d{3}_(\d{4}_\d{3}_)(\d{1,2})_(\d{1,4})(_*\..+)'
regex = re.compile(ooiregex)


def sort_function(filename):
    fname = os.path.basename(filename)
    match = regex.search(fname)
    if not match:
        logging.warning('File does not match dba file name')
        return
    else:
        mission_num = '{:02d}'.format(int(match.group(2)))
        segment_num = '{:04d}'.format(int(match.group(3)))
        sort_str = match.group(1) + mission_num + '_' + segment_num
    return sort_str
