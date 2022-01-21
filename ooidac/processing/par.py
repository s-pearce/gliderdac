from copy import deepcopy

import numpy as np
from ooidac.processing import logger, _add_nan_variable


def recalc_par(dba, sensor_dark, scale_factor):
    """ Recalculates PAR from the raw signal

    :param dba: GliderData instance
    :param sensor_dark: Dark Offset calibration parameter from the calibration
        certificate
    :param scale_factor: Appropriate units Wet Scale Factor calibration
        parameter from the calibration certificate
    :return: The GliderData instance with the new parameter added
    """
    if 'sci_bsipar_sensor_volts' not in dba.sensor_names:
        logger.warning(
            "sci_bsipar_sensor_volts is not present to recalculate PAR. "
            "Filling with NaNs instead.")
        # This needs to return a variable called `corrected_par` or else the
        # code will not continue with the rest of the potentially good data, so
        # this should return `corrected_par` full of nans since it can not be
        # corrected
        return _add_nan_variable(dba, "corrected_par", "sci_bsipar_par")
    par_volts = dba.getdata('sci_bsipar_sensor_volts')
    # remove the initialization where sensor volts == 0.0
    par_volts[par_volts == 0.0] = np.nan
    par_units = deepcopy(dba['sci_bsipar_par'])
    new_par = (par_volts - sensor_dark) / scale_factor
    par_units['data'] = new_par
    par_units['attrs']['comment'] = (
        "PAR recalculated from signal using calibration parameters")
    par_units['sensor_name'] = "corrected_par"
    par_units['attrs']['source_sensor'] = "sci_bsipar_sensor_volts"
    dba.add_data(par_units)

    return dba