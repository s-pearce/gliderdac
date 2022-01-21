import numpy as np
from ooidac.utilities import fwd_fill


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

    dba.add_data(pitch)
    dba.add_data(roll)

    return dba