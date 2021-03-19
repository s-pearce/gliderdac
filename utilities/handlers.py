import re


def backscatter_handler(attrs, proc=None):
    beta = _check_for_param(['beta'], attrs, proc) or 'sci_flbbcd_bb_units'
    wlngth = _check_for_param(
        ['wlngth', 'wavelength', 'radiation_wavelength'], attrs, proc) or 700
    if isinstance(wlngth, str):
        wlngth = int(re.search(r'(\d+)', wlngth).group(1))
    # put the radiation wavelength back in to attrs
    attrs['radiation_wavelength'] = wlngth
    attrs['radiation_wavelength_units'] = "nm"
    theta = _check_for_param(['theta'], attrs, proc) or 124.0
    xfactor = _check_for_param(['xfactor'], attrs, proc) or 1.076
    proc_dict = {
        "module": "fluorometer",
        "function": "glider_eco_backscatter",
        "data_args": {
            "beta": beta,
            "times": "m_present_time",
            "temp": "sci_water_temp",
            "salt": "salinity"
        },
        "param_args": {
            "wlngth": wlngth,
            "theta": theta,
            "xfactor": xfactor
        }
    }
    return proc_dict


def corrected_chlor_handler(attrs, proc=None):
    darkcnts = _check_for_param(
        ['counts_dark', 'dark_counts', 'dark_offset'],
        attrs, proc) or "<enter cal Dark Counts>"
    sf = _check_for_param(
        ['scale_factor'], attrs, proc) or "<enter cal Scale Factor>"

    proc_dict = {
        "module": "fluorometer",
        "function": "flo_chla",
        "data_args": {
            "counts_output": "sci_flbbcd_chlor_sig"
        },
        "param_args": {
            "counts_dark": darkcnts,
            "scale_factor": sf
        }
    }
    return proc_dict


def oxygen_handler(attrs, proc=None):
    # ToDo: a look up table?
    sref = (
        _check_for_param(["sref"], attrs, proc) or
        "<enter the salinity setting used if > 0 else remove this param>")
    # remove the old calibration coefficients since they are no longer needed
    attrs.pop("cal_coefs", None)

    proc_dict = {
        "module": "oxygen",
        "function": "aligned_salinity_correction",
        "data_args": {
            "molar_doxy": "sci_oxy4_oxygen",
            "timestamps": "llat_time",
            "ctd_pres": "llat_pressure",
            "ctd_temp": "sci_water_temp",
            "ctd_sp": "sci_water_cond",
            "lat": "llat_latitude",
            "lon": "llat_longitude"
        },
        "param_args": {
            "sref": sref
        }
    }
    return proc_dict


def corrected_par_handler(attrs, proc=None):
    sdark = _check_for_param(
        ['sensor_dark'], attrs, proc) or "<enter cal Sensor Dark>"
    # get the correct sensor_dark and scale_offset from the sdefs
    sf = _check_for_param(
        ['scale_factor'], attrs, proc) or "<enter cal Scale Factor>"
    proc_dict = {
        "module": "par",
        "function": "par_from_volts",
        "data_args": {
            "par_volts": "sci_bsipar_sensor_volts"
        },
        "param_args": {
            "sensor_dark": sdark,
            "scale_factor": sf
        }
    }
    return proc_dict


def pitch_handler(attrs, proc=None):
    proc_dict = {
        "module": "attitude",
        "function": "radtodegfill",
        "data_args": {
            "var": "m_pitch",
            "times": "llat_time"
        },
        "param_args": {
        }
    }
    return proc_dict


def roll_handler(attrs, proc=None):
    proc_dict = {
        "module": "attitude",
        "function": "radtodegfill",
        "data_args": {
            "var": "m_roll",
            "times": "llat_time"
        },
        "param_args": {
        }
    }
    return proc_dict


def sci_bbfl2s_chlor_sig_handler(attrs, proc=None):
    offset = _check_for_param([
        'dark_counts', 'dark_offset', 'clean_water_offset', 'cwo'],
        attrs, proc) or "<enter cal Clean Water Offset>"
    sf = _check_for_param(
        ['scale_factor'], attrs, proc) or "<enter cal Scale Factor>"
    proc_dict = {
        "module": "fluorometer",
        "function": "flo_scale_and_offset",
        "data_args": {
            "counts_output": "sci_bbfl2s_chlor_sig"
        },
        "param_args": {
            "offset": offset,
            "scale_factor": sf
        }
    }
    return proc_dict


def sci_bbfl2s_cdom_sig_handler(attrs, proc=None):
    offset = _check_for_param(
        ['dark_counts', 'dark_offset', 'clean_water_offset', 'cwo'],
        attrs, proc) or "<enter cal Clean Water Offset>"
    sf = _check_for_param(
        ['scale_factor'], attrs, proc) or "<enter cal Scale Factor>"
    proc_dict = {
        "module": "fluorometer",
        "function": "flo_scale_and_offset",
        "data_args": {
            "counts_output": "sci_bbfl2s_cdom_sig"
        },
        "param_args": {
            "offset": offset,
            "scale_factor": sf
        }
    }
    return proc_dict


def sci_bbfl2s_bb_sig_handler(attrs, proc=None):
    offset = _check_for_param(
        ['offset', 'dark_counts', 'dark_offset'],
        attrs, proc) or "<enter Cal Dark Counts>"
    sf = _check_for_param(
        ['scale_factor'], attrs, proc) or "<enter Cal Scale Factor>"
    # wavelength.  default to 660 for the bbfl2s
    wlngth = _check_for_param(
        ['wlngth', 'wavelength', 'radiation_wavelength'],
        attrs, proc) or 660
    # if the number is retrieved from a string to make it an int
    if isinstance(wlngth, str):
        wlngth = int(re.search(r'(d\+)', wlngth).group(1))

    # 124 is the default angle for a 3 channel ECO puck
    theta = _check_for_param(['theta'], attrs, proc) or 124.0
    # 1.076 is the default chi-factor for a 3 channel ECO puck
    xfactor = _check_for_param(['xfactor'], attrs, proc) or 1.076
    proc_dict = {
        "module": "fluorometer",
        "function": "glider_eco_backscatter_from_counts",
        "data_args": {
            "counts": "sci_bbfl2s_bb_sig",
            "times": "llat_time",
            "temp": "sci_water_temp",
            "salt": "salinity"
        },
        "param_args": {
            "offset": offset,
            "scale_factor": sf,
            "wlngth": wlngth,
            "theta": theta,
            "xfactor": xfactor
        }
    }
    return proc_dict


def rad_wlength_handler(attrs, proc=None):
    return None


processing_maps = {
    "backscatter": backscatter_handler,

    "corrected_chlor": corrected_chlor_handler,

    # get sref from ??  maybe a lookup table generated from the DAC spreadsheet
    "oxygen": oxygen_handler,
    "corrected_oxygen": oxygen_handler,
    "corrected_par": corrected_par_handler,
    "pitch": pitch_handler,
    "roll": roll_handler,
    "sci_bbfl2s_chlor_sig": sci_bbfl2s_chlor_sig_handler,
    "sci_bbfl2s_cdom_sig": sci_bbfl2s_cdom_sig_handler,
    "sci_bbfl2s_bb_sig": sci_bbfl2s_bb_sig_handler,
    "radiation_wavelength": rad_wlength_handler,
    "blank": {
        "comment":
            "fill out these fields and remove this comment for processing",
        "module": "",
        "function": "",
        "data_args": {},
        "param_args": {}
    }
}


def _check_for_param(paramnamelist, attrs, proc=None):
    """internal function to cycle through the possible choices for finding a
    parameter in the attrs or proc sensor_def dictionary"""
    param = None
    for paramname in paramnamelist:
        if paramname in attrs:
            param = attrs.pop(paramname)
        elif not param and proc and paramname in proc:
            param = proc.pop(paramname)
    return param
