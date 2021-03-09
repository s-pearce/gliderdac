import numpy as np
from ooidac.processing.general import fwd_fill


def pitch_and_roll(dba, fill='fwd fill'):
    """adds new sensors `pitch` and `roll` to a GliderData instance from
    `m_pitch` and `m_roll` converted to degrees from radians.

    :param dba:  A GliderData or DbaData instance
    :param fill:
    :return: dba:  The same GliderData instance with `pitch` and `roll` added
    """
    if fill == 'fwd fill':
        fill_function = fwd_fill
    elif fill == 'interp':
        def fill_function(param):
            filled_param = np.interp(
                dba.ts, dba.ts[np.isfinite(param)], param[np.isfinite(param)])
            return filled_param
    else:
        def fill_function(param):
            return param

    pitch = dba['m_pitch']
    pitch['sensor_name'] = 'pitch'
    pitch['attrs']['units'] = 'degrees'
    pitch['attrs']['comment'] = (
        'm_pitch converted to degrees and forward filled')
    pitch['data'] = fill_function(np.degrees(pitch['data']))

    roll = dba['m_roll']
    roll['sensor_name'] = 'roll'
    roll['attrs']['units'] = 'degrees'
    roll['attrs']['comment'] = (
        'm_roll converted to degrees and forward filled')
    roll['data'] = fwd_fill(np.degrees(roll['data']))

    dba.add_data(pitch['sensor_name'], pitch['data'],
                 **pitch['attrs'])
    dba.add_data(roll['sensor_name'], roll['data'],
                 **pitch['attrs'])

    return dba


def radtodegfill(var, times=None, fill='interp'):
    """Convert array from radians to degrees and fill NaNs

    Converts any of a Slocum glider's "attitude" variables `m_pitch`,
    `m_roll`, or `m_heading` from radians to degrees and fills to all of the
    timestamps with the chosen `fill` method

    :param var: array
        the attitude_rev variable to convert in radians
    :param times: array, only required for `fill`='interp'
        the timestamps or x variables to use with interpolation spacing.
        This is ignored if fill is not 'interp'.
    :param fill: ['fwd fill', 'interp', 'none'] default 'interp'
        The fill type used to fill NaN values.
    :return: array the same length as `attvar`
        the attitude_rev variables in degrees and filled to all
        timestamps
    """
    if fill == 'fwd fill':
        fill_function = fwd_fill
    elif fill == 'interp':
        assert isinstance(var, (list, tuple, np.ndarray)), (
            "input variable must be an array of values")
        if times is None:
            times = np.arange(len(var))
        assert len(var) == len(times), (
            "the `times` input is not the same length as `var`")

        def fill_function(param):
            filled_param = np.interp(
                times, times[np.isfinite(param)], param[np.isfinite(param)])
            return filled_param
    else:
        def fill_function(param):
            return param

    new_var = fill_function(np.degrees(var))

    return new_var
