import json
import sys
import os
from copy import deepcopy

from uts.handlers import processing_maps
from ooidac.readers.json_config import load


def convert_sdefs(file, nc_swap=True, file_dir=None):
    # sensor_defs.json rewrite
    # read json file into dictionary
    if not file_dir:
        file_dir = os.path.dirname(file)
    sdefs0 = load(file)
    sdef0keys = set(sdefs0.keys())
    defaults = {
        'llat_time', 'llat_pressure', 'llat_depth', 'llat_latitude',
        'llat_longitude', 'crs', 'platform', 'profile_id', 'profile_time',
        'profile_lat', 'profile_lon', 'sci_water_cond', 'sci_water_temp',
        'salinity', 'density', 'time_uv', 'lat_uv',
        'lon_uv', 'u', 'v'}

    ordered_sensor_defaults = [
        'platform', 'profile_id', 'profile_time', 'profile_lat', 'profile_lon',
        'time_uv', 'lat_uv', 'lon_uv', 'u', 'v', 'llat_time', 'llat_pressure',
        'llat_depth', 'llat_latitude', 'llat_longitude',
        'sci_water_cond', 'sci_water_temp', 'salinity', 'density']

    nondefaults = sdef0keys.difference(ordered_sensor_defaults)
    nondefaults = list(nondefaults)
    nondefaults.sort()

    # with open("resources/all_mdata_sensors.json", 'r') as fid:
    with open("../resources/all_mdata_sensors.json", 'r') as fid:
        mdata = json.load(fid)

    # reorder the entries starting with the llat_ variables, then the scalar
    # defaults, then the CTD defaults.

    # first reorder attributes in each sensor definition and add processing
    # sections for non defaults and known processing particles
    for key in sdefs0:
        sdef = sdefs0[key]
        dim = sdef.pop("dimension")
        isdim = sdef.pop("is_dimension", None)
        dimlen = sdef.pop("dimension_length", None)
        ncname = sdef.pop("nc_var_name")
        stype = sdef.pop("type", None)
        attrs = sdef.pop("attrs")
        proc = sdef.pop("processing", None)

        sdef["nc_var_name"] = ncname
        sdef["dimension"] = dim
        if stype:
            sdef["type"] = stype
        if isdim:
            sdef["is_dimension"] = isdim
        if dimlen:
            sdef["dimension_length"] = dimlen
        sdef["attrs"] = attrs

        # add processing section to non-default stuff using `processing_maps`
        # handler functions
        if proc and key not in processing_maps:
            updated_proc = processing_maps["blank"]
            updated_proc.update(proc)
            sdef["processing"] = updated_proc
        elif key in processing_maps:
            proc_handler = processing_maps[key]
            proc = proc_handler(attrs, proc)
            if proc:
                sdef["processing"] = proc
        elif (key not in defaults) and (key not in mdata):
            sdef["processing"] = processing_maps["blank"]

    # reorder the dict in a default way to reorder the entries.  Should I leave
    # out the defaults (at least the scalar ones?)
    sdefs1 = {}
    ordered_sensors = ordered_sensor_defaults
    ordered_sensors.extend(nondefaults)
    for sensor in ordered_sensors:
        sdefs1[sensor] = sdefs0[sensor]

    # write dict to a new json file.
    basefn = os.path.basename(file)
    basefn = os.path.splitext(basefn)[0]
    file2 = basefn + "_v2.json"
    # file2 = "sensor_defs_v2conv.json"
    with open(os.path.join(file_dir, file2), 'w') as fid:
        json.dump(sdefs1, fid, indent=4)

    # for future swap main dict key name to the NC variable name
    if nc_swap:
        sdefs2 = {}
        for key in sdefs1:
            newkey = sdefs1[key]['nc_var_name']
            sdefs2[newkey] = deepcopy(sdefs1[key])
            sdefs2[newkey]['source_sensor'] = key
            sdefs2[newkey].pop('nc_var_name')

        # ToDo: need a translation of llat_ and other variables in the
        #  processing list to their respective NC variable name (e.g.
        #  llat_latitude -> lat) and a second translation of intermediate
        #  variables back to their source_sensor
        #  (e.g. corrected_oxygen -> sci_oxy4_oxygen)

        # ToDo: figure out how for this version of definitions to pull from a
        #  particle that does not have a processing particle.  I'm thinking
        #  of either having an attribute that says "source_sensor" and the
        #  code will look for either "processing" or "source_sensor",
        #  or I force a processing particle for a non-default and instead of
        #  "module", "function", etc. only have "source_sensor" or if a
        #  singular value to keep too, I could make "value" a third option (
        #  e.g. for a radiation_wavelength type variable).

        with open(os.path.join(file_dir, 'v3_var_defs.json'), 'w') as fid:
            json.dump(sdefs2, fid, indent=4)


def my_run():
    file = "C:/Users/spearce/data/dac_tests/ce_326_R3/config/sensor_defs.json"
    file_dir = "C:/Users/spearce/temp/0_dac_json_write"
    convert_sdefs(file)


def main():
    file = sys.argv[1]
    file_dir = os.path.dirname(file)


if __name__ == "__main__":
    my_run()
