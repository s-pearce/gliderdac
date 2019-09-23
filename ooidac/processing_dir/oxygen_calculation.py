"""oxygen_calculation
Module of calculations and utilities for calculation of oxygen concentration
for an Aanderaa 4831 Oxygen Optode.

References
----------
AADI(2017). TD 269 OPERATING MANUAL: OXYGEN OPTODE 4330, 4831, 4835.
    https://www.aanderaa.com/media/pdfs/oxygen-optode-4330-4835-and-4831.pdf
Garcia, H.E. and Gordon, L.I. (1992). "Oxygen solubility in seawater:
    Better fitting equations". Limnol. Oceanogr. 37(6) 1307-1312.
Olson, Svein Rune. Aanderaa personal communication (2019)

Revision History
----------------
    2019-02-27: Stuart Pearce. Initial Code.
    2019-03-06: Stuart Pearce. Added ConcCoef adjustment parameters to O2 funcs
    2019-03-07: Stuart Pearce. Broke out the calculation steps into individual
        functions to use separately, and rewrote the main function in terms of
        the sub-functions.
"""
import numpy as np
import gsw

# Global definitions
KELVIN_OFFSET = 273.15  # The offset to convert temperature Celsius to Kelvin
ST_K = 298.15  # Standard Temperature in Kelvin. 25.0 deg C = 298.15 deg K


def calc_o2(cph, tmp, C, M, N, cc=(0, 1), S=0.0,
            NomAirPress=1013.25, NomAirMix=0.20946, idealgc=False):
    """Calculates oxygen concentration (uncorrected for salinity* or pressure)
    from raw data parameters calphase and temperature from an Aanderaa 4831
    optode using the MkII calculations for an optode that does not have the
    SVU (Stern-Volmer-Uchida) equation feature enabled.

    * if the Salinity calibration parameter within the optode is set to
    non-zero (e.g. 35), the MkII calculation partially corrects for the non-
    zero salinity as an average salinity.  The salinity and pressure correction
    algorithm that uses co-located CTD salinity still properly corrects for
    salinity if this is taken into account.

    Usage
    -----
    o2conc, airsat = calc_o2(
            calphase, temp, C, M, N, cc=[0 1],
            S=0.0, NomAirPress=1013.25, NomAirMix=0.20946 idealgc=False)

    Parameters
    ----------
    calphase : Calibrated_phase, Array/Scalar. [deg]
        Optode raw data output parameter CALPHASE
    temp : Temperature, Array/Scalar. [deg C]
        Optode thermistor TEMPERATURE measured near the foil
    C : Sensing Foil Coefficients. Array.
        Optode calibration parameters FOILCOEFA and FOILCOEFB concatenated
    M : Partial Pressure Polynomial Temperature Exponents. Array.
        Optode calibration parameters FOILPOLYDEGT
    N : Partial Pressure Polynomial Phase Exponents. Array.
        Optode calibration parameters FOILPOLYDEGO
    cc : Concentration adjustment coefficients. List.
        Optode calibration parameter CONCCOEF.  Default is [0, 1] if no
        CONCCOEF parameter is given (not included in early firmwares)
    S : Salinity configuration parameter. Scalar, Optional. [PSU]
        Default is 0.0. O2 concentration output from the opdtode uses the
        parameter SALINITY for this value in the calculation.
        If it is desired to calculate the oxygen exactly the same as the
        optode did and the SALINITY parameter was non-zero, set this.
    NomAirPress : Nominal Air Pressure, Scalar, Optional.  [hPa]
        Default is 1013.25. O2 concentration output from the optode uses
        the parameter NOMAIRPRESS for this value in the calculation.
    NomAirMix : Nominal Air Mixture, Scalar, Optional.
        The mole fraction of O2 in dry air.  Default is 0.20946.
        O2 concentration output from the optode uses the parameter
        NOMAIRMIX for this value in the calculation.
    idealgc : Bool, optional.
        Use the ideal gas constant (44.615) instead of the real (44.659).
        True will use the ideal gas constant for backward compatibility with
        optode firmwares 1.*.* before 2012. Default is False (use the real gas
        constant)

    Returns
    -------
    A tuple of o2conc and airsat.

    o2conc : Oxygen concentration, Array/Scalar. [micro-moles/L]
        O2 concentration uncorrected for salinity and pressure.
    airsat : Air saturation, Array/Scalar. [%]

    References
    ----------
        AADI(2017). TD 269 OPERATING MANUAL: OXYGEN OPTODE 4330, 4831, 4835.
            https://www.aanderaa.com/media/pdfs/
                oxygen-optode-4330-4835-and-4831.pdf
        Garcia, H.E. and Gordon, L.I. (1992). "Oxygen solubility in seawater:
            Better fitting equations". Limnol. Oceanogr. 37(6) 1307-1312.

    Revision History
    ----------------
        2019-02-27: Stuart Pearce. Initial Code.
        2019-03-06: Stuart Pearce. Added ConcCoef adjustment parameter inputs
    """
    # O2 partial pressure
    partialpres = partial_pressure(cph, tmp, C, M, N)

    # vapour pressure pvapor(t)
    pvapor = vapor_pressure(tmp)

    # air saturation percentage
    airsaturation = air_saturation(
            partialpres, pvapor, NomAirPress=NomAirPress, NomAirMix=NomAirMix)

    if idealgc:
        gas_const = 44.615
    else:
        gas_const = 44.659

    # Oxygen solubility (cm3/dm3) from Garcia and Gordon (1992)
    oxysol = oxygen_solubility(tmp, S=S) * gas_const

    # O2 concentration calculated as micro-mole/L
    o2Conc = oxysol * airsaturation / 100.

    # Calibration adjustment using the ConcCoef calibration coefficents
    o2Conc = cc[0] + o2Conc*cc[1]

    return o2Conc, airsaturation


def calphase_from_rph(c1rph, c2rph, phaseCoefs):
    """
    Recalculates calphase for an Aanderaa 4831 oxygen optode from the raw blue
    and red excitation light phases using the calibration parameters determined
    for the optode.

    Parameters
    ----------
    c1rph : Blue excitation light phase, Array/Scalar. [deg]
            Optode raw data output parameter C1RPH
    c2rph : Red excitation light phase, Array/Scalar. [deg]
            Optode raw data output parameter C2RPH
    phaseCoefs : Coefficients for calculating CALPHASE from TCPHASE. List.
            Optode calibration parameter PHASECOEFS

    Returns
    -------
    calphase : Calibrated_phase, Array/Scalar. [deg]
            A recalculation of Optode raw data output parameter CALPHASE

    References
    ----------
    AADI(2017). TD 269 OPERATING MANUAL: OXYGEN OPTODE 4330, 4831, 4835.
        https://www.aanderaa.com/media/pdfs/
            oxygen-optode-4330-4835-and-4831.pdf
    """
    tcphase = c1rph - c2rph
    calphase = (
            phaseCoefs[0]
            + phaseCoefs[1] * tcphase
            + phaseCoefs[2] * tcphase**2
            + phaseCoefs[3] * tcphase**3
            )
    return calphase


def temp_from_rawtemp(rawtemp, tempCoefs):
    """
    Recalculates temperature for an Aanderaa 4831 oxygen optode from the raw
    temperature thermistor voltage using the calibration parameters determined
    for the optode.

    Parameters
    ----------
    rawtemp : Raw temperature themistor reading. Array/Scalar. [mV]
            Optode raw data output parameter RAWTEMP
    tempCoefs : Coefficients for calculating TEMPERATURE from RAWTEMP. List.
            Optode calibration parameters TEMPCOEFS found in the calibration
            certificate but not in the "get all" printout.

    Returns
    -------
    temp : Temperature, Array/Scalar. [deg C]
            A recalculation of Optode thermistor TEMPERATURE measured near
            the foil

    References
    ----------
    AADI(2017). TD 269 OPERATING MANUAL: OXYGEN OPTODE 4330, 4831, 4835.
        https://www.aanderaa.com/media/pdfs/
            oxygen-optode-4330-4835-and-4831.pdf
    """
    temp = (
            tempCoefs[0]
            + tempCoefs[1] * rawtemp
            + tempCoefs[2] * rawtemp**2
            + tempCoefs[3] * rawtemp**3
            )
    return temp


def partial_pressure(cph, tmp, C, M, N):
    """
    Calculates partial pressure [hPa] from calphase and temperature for an
    Aanderaa 4831 oxygen optode using calibration parameters determined for the
    optode.

    Usage
    -----
    PP = partial_pressure(calphase, temp, C, M, N)

    Parameters
    ----------
    calphase = array-like or scalar. [deg]
        CalPhase, from Optode raw data output
    temp = array-like or scalar. [deg C]
        Temperature, from Optode data output
    C = array-like.
        Optode calibration coefficients *FoilCoefA* and `FoilCoefB`
        concatenated.
    M = array-like
        Optode calibration coefficient *FoilPolyDegT*
    N = array-like.
        Optode calibration coefficient `FoilPolyDegO`

    Returns
    -------
    PP = array-like or scalar. [hPa]
            Partial pressure.

    References
    ----------
        AADI(2017). TD 269 OPERATING MANUAL: OXYGEN OPTODE 4330, 4831, 4835.
            https://www.aanderaa.com/media/pdfs/
                oxygen-optode-4330-4835-and-4831.pdf

    Revision History
    ----------------
        2019-03-07: Stuart Pearce. Initial Code.
    """
    C = np.atleast_1d(C)
    M = np.atleast_1d(M)
    N = np.atleast_1d(N)
    # initialize container to be the same size as calphase
    partialpress = np.zeros_like(cph)
    # Partial pressure is a polynomial summation calculated here in a for loop
    for ii in range(len(C)):
        partialpress = partialpress + C[ii] * tmp**M[ii] * cph**N[ii]
    return partialpress


def vapor_pressure(temp):
    """
    Calculates vapor pressure [hPa] from temperature.

    Parameters
    ----------
    temp = array-like or scalar. [deg C]
        Temperature, from Optode data output

    Returns
    -------
    vpress = array-like or scalar. [hPa]
            vapor pressure.

    References
    ----------
        AADI(2017). TD 269 OPERATING MANUAL: OXYGEN OPTODE 4330, 4831, 4835.
            https://www.aanderaa.com/media/pdfs/
                oxygen-optode-4330-4835-and-4831.pdf
    """
    temp = np.atleast_1d(temp)
    # KELVIN_OFFSET is defined in the Global definitions section at the top
    # of the module.

    # vapour pressure pvapor(t)
    pvapor = np.exp(
            52.57
            - (6690.9/(temp + KELVIN_OFFSET))
            - 4.6810 * np.log(temp + KELVIN_OFFSET))
    return pvapor


def air_saturation(partialpress, pvapor,
                   NomAirPress=1013.25, NomAirMix=0.20946):
    """Calculates air saturation for an Aanderaa 4831 oxygen optode from
    calculated partial pressure [hPa] and vapor pressure [hPa]

    Parameters
    ----------
    partialpress : array-like or scalar. [hPa]
        Partial pressure, calculated from Optode CalPhase and Temperature.
    pvapor : array-like or scalar. [hPa]
        Vapor Pressure, calculated from Optode Temperature.
    NomAirPres : scalar, optional.  [hPa]
        Nominal Air Pressure.  Optode configuration parameter *NomAirPress*,
        which is used in the calculation of air saturation (and thus O2
        concentration) in the optode. Default is 1013.25.
    NomAirMix = scalar, optional.
        Nominal Air Mixture.  The mole fraction of O2 in dry air. Optode
        configuration parameter *NomAirMix*, which is used in the calculation
        of air saturation (and thus O2 concentration) in the optode. Default is
        0.20946.

    Returns
    -------
    airsat = array-like or scalar. [%]
        Air saturation.
    """
    partialpress = np.atleast_1d(partialpress)
    pvapor = np.atleast_1d(pvapor)
    # air saturation percent
    airsaturation = (partialpress * 100.) / ((NomAirPress-pvapor) * NomAirMix)
    return airsaturation


def oxygen_solubility(temp, S=0.0):
    """Calculates oxygen solubility from temperature and salinity for
    calculation of O2 concentration from an oxygen sensor using the Garcia and
    Gordon equation (1992).

    Parameters
    ----------
    temp : array-like or scalar. [deg C]
        Temperature, from Optode data output or a co-located CTD.
    S :  array-like or scalar, optional. [PSU]
        Salinity, from the Optode configuration parameter or a co-located CTD.
        O2 concentration output from the optode uses the configuration
        parameter *Salinity* in its calculation of oxygen solubility. To
        replicate the calculation within the optode, set this to the same value
        as the configuration parameter.  Otherwise, exact salinities from a
        co-located CTD can be used.
        Default is 0.0.

    Returns
    -------
    oxysol = array-like or scalar. [cm3/dm3]
        Oxygen Solubility.

    Reference
    ---------
    Garcia, H.E. and Gordon, L.I. (1992). "Oxygen solubility in seawater:
        Better fitting equations". Limnol. Oceanogr. 37(6) 1307-1312.
    """
    temp = np.atleast_1d(temp)
    # KELVIN_OFFSET and ST_K (Standard Temperature in Kelvin) are defined in
    # the Global definitions section at the top of the module.

    # scaled temperature
    Ts = np.log((ST_K - temp) / (KELVIN_OFFSET + temp))

    # Garcia and Gordon coefficients
    A = [2.00856, 3.22400, 3.99063, 4.80299, 9.78188e-1, 1.71069]
    B = [-6.24097e-3, -6.93498e-3, -6.90358e-3, -4.29155e-3]
    C0 = -3.11680e-7

    oxysol = np.exp(  # Oxygen solubility (cm3/dm3)
            A[0] + A[1]*Ts + A[2]*Ts**2 + A[3]*Ts**3 + A[4]*Ts**4 + A[5]*Ts**5
            + S * (B[0] + B[1]*Ts + B[2]*Ts**2 + B[3]*Ts**3)
            + C0*S**2)
    return oxysol


def do2_SVU(calphase, temp, csv, conc_coef=np.array([0.0, 1.0]), salt=0):
    """
    Calculates oxygen concentration (uncorrected for salinity* or pressure)
    from raw data parameters calphase and temperature from an Aanderaa 4831
    optode using the SVU calculations for an optode that does have the
    SVU (Stern-Volmer-Uchida) equation feature enabled.

    * if the Salinity calibration parameter within the optode is set to
    non-zero (e.g. 35), the SVU calculation partially corrects for the non-
    zero salinity as an average salinity.  The salinity and pressure correction
    algorithm that uses co-located CTD salinity still properly corrects for
    salinity if this is taken into account.

    Usage
    -----
    DO = do2_SVU(calphase, temp, csv, conc_coef, salt)

    Parameters
    ----------
    calphase : calibrated phase from an Oxygen sensor [deg], DOCONCS-DEG_L0
        (see DOCONCS DPS)
    temp : oxygen sensor foil temperature T(optode) [deg C], (see DOCONCS DPS)
    csv : Stern-Volmer-Uchida Calibration Coefficients array.
        7 element float array, (see DOCONCS DPS)
    conc_coef : 'secondary' calibration coefficients: an array of offset and
        slope coefficients to apply to the result of the SVU equation.
        See Notes.
        conc_coef[0, 0] = offset
        conc_coef[0, 1] = slope
    salt : preset salinity parameter on the optode

    Returns
    -------
    o2conc : Oxygen concentration, Array/Scalar. [micro-moles/L]
        dissolved oxygen [micro-mole/L], DOCONCS_L1. see Notes.
        O2 concentration uncorrected for salinity and pressure.

    References
    ----------
        AADI(2017). TD 269 OPERATING MANUAL: OXYGEN OPTODE 4330, 4831, 4835.
            https://www.aanderaa.com/media/pdfs/
                oxygen-optode-4330-4835-and-4831.pdf
        Garcia, H.E. and Gordon, L.I. (1992). "Oxygen solubility in seawater:
            Better fitting equations". Limnol. Oceanogr. 37(6) 1307-1312.

    Revision History
    ----------------
        2019-02-27: Stuart Pearce. Initial Code.
        2019-03-06: Stuart Pearce. Added ConcCoef adjustment parameter inputs


    Usage:

        DO = do2_SVU(calphase, temp, csv, conc_coef)

            where



    Example:
        csv = np.array([0.002848, 0.000114, 1.51e-6, 70.42301, -0.10302,
                        -12.9462, 1.265377])
        calphase = 27.799
        temp = 19.841

        DO = do2_SVU(calphase, temp, csv)
        print DO
        > 363.900534505

    Implemented by:
        2013-04-26: Stuart Pearce. Initial Code.
        2015-04-10: Russell Desiderio. Revised code to work with CI implementation
                    of calibration coefficients: they are to be implemented as time-
                    vectorized arguments (tiled in the time dimension to match the
                    number of data packets). Fix for "blocker #2972".
        2015-08-04: Russell Desiderio. Added documentation.
        2015-08-10: Russell Desiderio. Added conc_coef calibration array to argument list.
                    Required to be a 2D row vector for broadcasting purposes.
        2015-10-28: Russell Desiderio. Added conc_coef = np.atleast_2d(conc_coef) line so
                    that function will now accept conc_coef as a 1D array (so that 1D array
                    entries in Omaha cal sheets won't result in DPA exceptions being raised).
                    So. Also changed default value for conc_coef in argument list to be
                    the 1D array [0.0, 1.0].

    Notes:

        General:

            The DOCONCS_L1 data product has units of micromole/liter; SAF incorrectly
            lists the units for this L1 product as micromole/kg. (To change units from
            mmole/L to mmole/kg, salinity is required, making the result an L2 data
            product).

            The DOCONCS_L1 data product is uncorrected for salinity and pressure.

        Temperature dependence:

            The optode sensor's thermistor temperature should be used whenever possible
            because for the OOI DOSTAs (model 4831) it is situated directly at the sensor
            foil and the SVU cal coefficients are derived in part to compensate for the
            change in oxygen permeability through the foil as a function of its temperature.

            The time constant of the model 4831 thermistor is < 2 seconds. Because the foil
            and therefore the calphase response time itself is 8 sec or 24 sec depending on
            the particular optode, there is little or no advantage to be gained by using a
            temperature sensor (eg, from a CTD) with a faster response. It is better to make
            sure that the temperature used most accurately reflects the foil temperature.

            On gliders, there is often a difference in CTD and optode temperature readings of
            1 degree Celsius, which translates to about a 5% difference in calculated oxygen
            concentration for a range of typical water column conditions.

        Conc_coef (this information is not currently in the DPS):

            Aanderaa uses two calibration procedures for the 4831 optode. The primary 'multi-point'
            calibration, done in Norway, determines the SVU foil coefficients (variable csv in the
            DPA). The secondary two-point calibration, done in Ohio, corrects the multi-point
            calibration calculation against 0% oxygen and 100% oxygen data points to provide the
            conc_coef values. (Aanderaa is in the process of changing the secondary cal to a one
            point cal, using just the 100% oxygen data point, but the result will still be expressed
            as offset and slope conc_coef values.) For standard optode refurbishment Aanderaa recommends
            a secondary calibration instead of a new multi-point SVU foil calibration.

            Secondary calibrations are not done on new optodes nor on optodes with new determinations
            of the SVU foil coefficients; in these cases Aanderaa sets the conc_coef values to 0 (offset)
            and 1 (slope) in the optode firmware by default. Conc_coef determinations resulting from the
            secondary calibration procedure are also incorporated into the optode firmware and are also
            listed on the Aanderaa Form No. 710 calibration certificate, although they are currently
            mislabelled on this form as "PhaseCoef".

            The conc_coef correction to optode-calculated values for oxygen concentration is automatically
            applied by the optode firmware. However, this correction must be done manually when oxygen
            concentration is calculated from calphase and optode temperature external to the optode, as in
            this DPA do2_SVU.

    References:

        OOI (2012). Data Product Specification for Oxygen Concentration
            from "Stable" Instruments. Document Control Number
            1341-00520. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level
            >> 1341-00520_Data_Product_SPEC_DOCONCS_OOI.pdf)

        Aanderaa Data Instruments (August 2012). TD 269 Operating Manual Oxygen Optode 4330, 4831, 4835.

        August 2015. Shawn Sneddon, Xylem-Aanderaa technical support, MA, USA, 800-765-4974
    """
    conc_coef = np.atleast_2d(conc_coef)
    # this will work for both old and new CI implementations of cal coeffs.
    csv = np.atleast_2d(csv)

    # Calculate DO using Stern-Volmer:
    Ksv = csv[:, 0] + csv[:, 1]*temp + csv[:, 2]*(temp**2)
    P0 = csv[:, 3] + csv[:, 4]*temp
    Pc = csv[:, 5] + csv[:, 6]*calphase
    DO = ((P0/Pc) - 1) / Ksv

    # scaled temperature
    Ts = np.log((ST_K - temp) / (KELVIN_OFFSET + temp))

    # Garcia and Gordon coefficients
    B = [-6.24097e-3, -6.93498e-3, -6.90358e-3, -4.29155e-3]
    C0 = -3.11680e-7

    oxysol = np.exp(  # Oxygen solubility (cm3/dm3)
            + salt * (B[0] + B[1]*Ts + B[2]*Ts**2 + B[3]*Ts**3)
            + C0*salt**2)

    DO = DO * oxysol

    # apply refurbishment calibration
    # conc_coef can be a 2D array of either 1 row or DO.size rows.
    DO = conc_coef[:, 0] + conc_coef[:, 1] * DO
    return DO


def do2_salinity_correction(DO, P, T, SP, lat, lon, pref=0):
    """
    Description:

        Calculates the data product DOXYGEN_L2 (renamed from DOCONCS_L2) from DOSTA
        (Aanderaa) instruments by correcting the the DOCONCS_L1 data product for
        salinity and pressure effects and changing units.

    Usage:

        DOc = do2_salinity_correction(DO,P,T,SP,lat,lon, pref=0)

            where

        DOc = corrected dissolved oxygen [micro-mole/kg], DOXYGEN_L2
        DO = uncorrected dissolved oxygen [micro-mole/L], DOCONCS_L1
        P = PRESWAT water pressure [dbar]. (see
            1341-00020_Data_Product_Spec_PRESWAT). Interpolated to the
            same timestamp as DO.
        T = TEMPWAT water temperature [deg C]. (see
            1341-00010_Data_Product_Spec_TEMPWAT). Interpolated to the
            same timestamp as DO.
        SP = PRACSAL practical salinity [unitless]. (see
            1341-00040_Data_Product_Spec_PRACSAL)
        lat, lon = latitude and longitude of the instrument [degrees].
        pref = pressure reference level for potential density [dbar].
            The default is 0 dbar.

    Example:
        DO = 433.88488978325478
        do_t = 1.97
        P = 5.4000000000000004
        T = 1.97
        SP = 33.716000000000001
        lat,lon = -52.82, 87.64

        DOc = do2_salinity_correction(DO,P,T,SP,lat,lon, pref=0)
        print DO
        > 335.967894709

    Implemented by:
        2013-04-26: Stuart Pearce. Initial Code.
        2015-08-04: Russell Desiderio. Added Garcia-Gordon reference.

    References:
        OOI (2012). Data Product Specification for Oxygen Concentration
            from "Stable" Instruments. Document Control Number
            1341-00520. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level
            >> 1341-00520_Data_Product_SPEC_DOCONCS_OOI.pdf)

        "Oxygen solubility in seawater: Better fitting equations", 1992,
        Garcia, H.E. and Gordon, L.I. Limnol. Oceanogr. 37(6) 1307-1312.
        Table 1, 5th column.
    """

    # density calculation from GSW toolbox
    SA = gsw.SA_from_SP(SP, P, lon, lat)
    CT = gsw.CT_from_t(SA, T, P)
    pdens = gsw.rho(SA, CT, pref)  # potential referenced to p=0

    # Convert from volume to mass units:
    DO = 1000*DO/pdens

    # Pressure correction:
    DO = (1 + (0.032*P)/1000) * DO

    # Salinity correction (Garcia and Gordon, 1992, combined fit):
    S0 = 0
    ts = np.log((298.15-T)/(273.15+T))
    B0 = -6.24097e-3
    B1 = -6.93498e-3
    B2 = -6.90358e-3
    B3 = -4.29155e-3
    C0 = -3.11680e-7
    Bts = B0 + B1*ts + B2*ts**2 + B3*ts**3
    DO = np.exp((SP-S0)*Bts + C0*(SP**2-S0**2)) * DO
    return DO
