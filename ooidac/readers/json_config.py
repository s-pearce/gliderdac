import json
"""json_config.py
A module for reading and writing the JSON configuration files associated with
Glider DAC deployment runs.
"""


def load(path):
    """JSON loader that allows for comments

    The JSON configuration files required for Glider DAC deployment processing
    can have JavaScript style line comments that begin with a //.  Comments
    must be on their own line and not in line with other JSON elements, but may
    have any amount of leading whitespace for indentation purposes.

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
    lineno = 0
    for line in lines:
        if line.strip().startswith("//"):
            comment_lines.append((lineno, line))
        lineno += 1

    jsonlines = lines.copy()

    for comment in comment_lines:
        jsonlines.remove(comment[1])

    # Error control here is to adjust JSONDecoderError to account for removed
    # comment lines
    try:
        json_dict = json.loads("".join(jsonlines))
    except json.JSONDecodeError as e:
        # adjust the line number and character numbers to count the removed
        # comment lines
        lineno = e.lineno
        num_chars = e.pos
        for commentlineno, comment in comment_lines:
            if commentlineno <= lineno:
                lineno += 1
                num_chars += len(comment)
            else:
                break

        # adjust the error message with new values
        newmsg = e.msg + (
            ': line {:d} column {:d} (char {:d}) in:\n  {:s}\n'.format(
                lineno, e.colno, num_chars, path)
        )
        # to print the actual lines, we must convert lineno to the lines list
        # index by subtracting 1 due to counting starting from 0 while
        # numbers in a text editor will start counting from 1
        lineno -= 1
        for ii in range(lineno-1, lineno+2):
            if ii == lineno:
                arrow = '->'
            else:
                arrow = '  '
            newmsg += '{:s}{:d}:{:s}'.format(arrow, ii+1, lines[ii])
        e.args = (newmsg, )
        raise
    return json_dict
