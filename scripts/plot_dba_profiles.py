#!/usr/bin/env python

import logging
import os
import argparse
import sys
import glob

codepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, codepath)

from ooidac.plotting import plot_profiles_for_dba


def main(args):
    """Parse the specified Slocum glider dba file and plot the yo and the
    indexed profiles"""

    # Set up logger
    log_level = getattr(logging, args.loglevel.upper())
    log_format = '%(module)s:%(levelname)s:%(message)s [line %(lineno)d]'
    logging.basicConfig(format=log_format, level=log_level)

    dba_files = glob.glob(args.dba_files)

    for dba_file in dba_files:

        plot_profiles_for_dba(dba_file, args.unfiltered)

    return 0


if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser(
        description=main.__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument(
        'dba_files',
        help='Path to a dba file or files')

    arg_parser.add_argument(
        '-u', '--unfiltered',
        help='Do not filter the indexed profiles',
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
