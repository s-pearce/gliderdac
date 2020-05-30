import os
import json
"""json_config.py
A module for reading and writing the JSON configuration files associated with
Glider DAC deployment runs.
"""


def load(path):
    """JSON loader that allows for comments

    The JSON configuration files required for Glider DAC deployment processing
    can have line comments that begin with a #.  Comments must be on their
    own line and not in line with other JSON elements, but may have any amount
    of leading whitespace for indentation purposes.

    :param path: str
        The string path to the JSON config file (e.g. sensor_defs.json,
        deployment.json, etc.)
    :return: dictionary
        The JSON object returned as a dictionary.
    """
    # No error control included here to allow the native FileNotFoundError
    # exception to raise if path is not valid.
    with open(path, 'r') as fid:
        lines = fid.readlines()
    comment_lines = []
    for line in lines:
        if line.strip().startswith("#"):
            comment_lines.append(line)

    for comment in comment_lines:
        lines.remove(comment)

    # No error control here either to allow native JSONDecoder exceptions
    return json.loads("".join(lines))