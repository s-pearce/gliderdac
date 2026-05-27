import numpy as np
from scipy import interpolate


def water_depth(prof):
    """Add `water_depth` and `altitude` to a GliderData 
    profile instance by foward filling `m_water_depth` and
    `m_altitude` to all timestamps

    :param prof:  A GliderData or DbaData profile instance
    :return: dba:  The same GliderData instance with `pitch` and `roll` added
    """

    wdepth = prof['m_water_depth']
    wdepth['sensor_name'] = 'water_depth'
    wdepth['attrs']['comment'] = (
        'm_water_depth forward filled to science timestamps')

    alt = prof['m_altitude']
    alt['sensor_name'] = 'altitude'
    alt['attrs']['comment'] = (
        'm_altitude forward filled to science timestamps')

    ## Forward fill water depth and altitude to all timestamps
    ## so that when the profile is reduced to valid science timestamps
    ## water depth and altitude are kept.

    # first set any -1 water depths to NaN for both water depth and
    # altimeter
    n = len(wdepth['data'])
    wd_neg = wdepth['data'] == -1
    wdepth['data'][wd_neg] = np.nan
    alt['data'][wd_neg] = np.nan
    
    # then find non-nan indices
    wd_ii = np.isfinite(wdepth['data'])
    alt_ii = np.isfinite(alt['data'])

    # forward-fill (same as "previous") to all profile timestamps
    # We only are interested in profiles where number of altimeter hits
    # are greater than 2.  Less than 2 is not reliable.
    if np.sum(wd_ii) > 2 and np.sum(alt_ii) > 2:
        
        sts = prof['sci_m_present_time']['data']
        sts_ii = np.isfinite(sts)
        
        last_wd_val = wdepth['data'][wd_ii][-1]
        last_alt_val = alt['data'][alt_ii][-1]
        
        tmp_ff_wd = interpolate.interp1d(
            prof.ts[wd_ii], wdepth['data'][wd_ii], kind='previous',
            bounds_error=False, fill_value=(np.nan, last_wd_val))(sts[sts_ii])
        ff_wd = np.full(n, np.nan)
        ff_wd[sts_ii] = tmp_ff_wd
        wdepth['data'] = ff_wd
        
        tmp_ff_alt = interpolate.interp1d(
            prof.ts[alt_ii], alt['data'][alt_ii], kind='previous',
            bounds_error=False, fill_value=(np.nan, last_alt_val))(sts[sts_ii])
        ff_alt = np.full(n, np.nan)
        ff_alt[sts_ii] = tmp_ff_alt
        alt['data'] = ff_alt
    else:
        wdepth['data'] = np.full(n, np.nan)
        alt['data'] = np.full(n, np.nan)

    prof.add_data(wdepth)
    prof.add_data(alt)

    return prof