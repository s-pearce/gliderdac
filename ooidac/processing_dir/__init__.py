from .fluorometer import flo_chla


def recalc_flbbcd(dba, sensor, dark_offset, scale_factor):
    """ Recalculates any of the flbbcd units from the raw signal given the
    calibration parameters from either the calibration certificate document
    or determined independently and inserts it back into the data set

    :param dba: GliderData instance
    :param sensor: the name of the flbbcd sensor to be re-calculated
    :param dark_offset: Dark Offset calibration parameter
    :param scale_factor: Scale Factor calibration parameter
    :return: The GliderData instance with the new parameter added
    """
    signal = dba.getdata(sensor)
    particle = dba[sensor]
    new_units = flo_chla(signal, dark_offset, scale_factor)
    particle['data'] = new_units
    particle['comment'] = (
        "Recalculated from raw signal")
    particle['sensor_name'] = "corrected_" + sensor
    dba.add_data(particle)

    return dba
