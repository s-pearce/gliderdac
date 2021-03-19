
from configuration import (
    PROC_PRES_VAR, PROC_DEPTH_VAR, PROC_TIME_VAR,
    PROC_LON_VAR, PROC_LAT_VAR)
nan = float('nan')

default_sdefs = {
    # rest should be filled out by the user in the configurations
    "platform": {
        "nc_var_name": "platform",
        "dimension": None,
        "type": "i4",
        "attrs": {
            "ancillary_variables": " ",
            "instrument": (
                "instrument_ctd")
        }
    },
    "profile_id": {
        "nc_var_name": "profile_id",
        "dimension": None,
        "type": "i4",
        "attrs": {
            "_FillValue": -999,
            "ancillary_variables": "profile_time",
            "comment": (
                "Sequential profile number within the trajectory. "
                "This value is unique in each file that is part of a single "
                "trajectory/deployment."),
            "long_name": "Profile ID",
            "valid_min": 1
        }
    },
    "profile_time": {
        "nc_var_name": "profile_time",
        "dimension": None,
        "type": "f8",
        "attrs": {
            "_FillValue": -999.0,
            "axis": "T",
            "calendar": "gregorian",
            "comment": (
                "Timestamp corresponding to the mid-point of the profile"),
            "long_name": "Profile Center Time",
            "observation_type": "calculated",
            "platform": "platform",
            "standard_name": "time",
            "units": "seconds since 1970-01-01T00:00:00Z"
        }
    },
    "profile_lat": {
        "nc_var_name": "profile_lat",
        "dimension": None,
        "type": "f8",
        "attrs": {
            "_FillValue": -999.0,
            "comment": (
                "Value is interpolated to provide an estimate of the "
                "latitude at the mid-point of the profile"),
            "units": "degrees_north",
            "platform": "platform",
            "coordinate_reference_frame": "urn:ogc:crs:EPSG::4326",
            "reference": "WGS84",
            "axis": "Y",
            "observation_type": "calculated",
            "long_name": "Profile Center Latitude",
            "valid_min": -90.0,
            "valid_max": 90.0,
            "standard_name": "latitude"
        }
    },
    "profile_lon": {
        "nc_var_name": "profile_lon",
        "dimension": None,
        "type": "f8",
        "attrs": {
            "_FillValue": -999.0,
            "comment": (
                "Value is interpolated to provide an estimate of the "
                "longitude at the mid-point of the profile"),
            "units": "degrees_east",
            "platform": "platform",
            "coordinate_reference_frame": "urn:ogc:crs:EPSG::4326",
            "reference": "WGS84",
            "axis": "X",
            "observation_type": "calculated",
            "long_name": "Profile Center Longitude",
            "valid_min": -180.0,
            "valid_max": 180.0,
            "standard_name": "longitude"
        }
    },
    "time_uv": {
        "nc_var_name": "time_uv",
        "dimension": None,
        "type": "f8",
        "attrs": {
            "_FillValue": nan,
            "long_name": "Depth-Averaged Time",
            "standard_name": "time",
            "calendar": "gregorian",
            "comment": (
                "The depth-averaged current is an estimate of the net current "
                "measured while the glider is underwater. The value is "
                "calculated over the entire underwater segment, which may "
                "consist of 1 or more dives."),
            "observation_type": "calculated",
            "units": "seconds since 1970-01-01T00:00:00Z"
        }
    },
    "lat_uv": {
        "nc_var_name": "lat_uv",
        "dimension": None,
        "type": "f8",
        "attrs": {
            "_FillValue": nan,
            "long_name": "Depth-Averaged Latitude",
            "standard_name": "latitude",
            "comment": (
                "The depth-averaged current is an estimate of the net current "
                "measured while the glider is underwater. The value is "
                "calculated over the entire underwater segment, which may "
                "consist of 1 or more dives."),
            "observation_type": "calculated",
            "platform": "platform",
            "units": "degrees_north",
            "valid_max": 90.0,
            "valid_min": -90.0
        }
    },
    "lon_uv": {
        "nc_var_name": "lon_uv",
        "dimension": None,
        "type": "f8",
        "attrs": {
            "_FillValue": nan,
            "long_name": "Depth-Averaged Longitude",
            "standard_name": "longitude",
            "comment": (
                "The depth-averaged current is an estimate of the net current "
                "measured while the glider is underwater. The value is "
                "calculated over the entire underwater segment, which may "
                "consist of 1 or more dives."),
            "observation_type": "calculated",
            "platform": "platform",
            "units": "degrees_east",
            "valid_max": 180.0,
            "valid_min": -180.0
        }
    },
    "u": {
        "nc_var_name": "u",
        "dimension": None,
        "type": "f8",
        "attrs": {
            "_FillValue": nan,
            "long_name": "Depth-Averaged Eastward Sea Water Velocity",
            "standard_name": "eastward_sea_water_velocity",
            "comment": (
                "The depth-averaged current is an estimate of the net current "
                "measured while the glider is underwater. The value is "
                "calculated over the entire underwater segment, which may "
                "consist of 1 or more dives."),
            "observation_type": "calculated",
            "platform": "platform",
            "units": "m s-1",
            "valid_max": 10.0,
            "valid_min": -10.0
        }
    },
    "v": {
        "nc_var_name": "v",
        "dimension": None,
        "type": "f8",
        "attrs": {
            "_FillValue": nan,
            "long_name": "Depth-Averaged Northward Sea Water Velocity",
            "standard_name": "northward_sea_water_velocity",
            "comment": (
                "The depth-averaged current is an estimate of the net current "
                "measured while the glider is underwater. The value is "
                "calculated over the entire underwater segment, which may "
                "consist of 1 or more dives."),
            "observation_type": "calculated",
            "platform": "platform",
            "units": "m s-1",
            "valid_max": 10.0,
            "valid_min": -10.0
        }
    },
    PROC_TIME_VAR: {
        "nc_var_name": "time",
        "dimension": "time",
        "type": "f8",
        "is_dimension": True,
        "attrs": {
            "axis": "T",
            "calendar": "gregorian",
            "standard_name": "time",
            "long_name": "Time",
            "units": "seconds since 1970-01-01T00:00:00Z",
            "_FillValue": nan,
            "ancillary_variables": " ",
            "observation_type": "measured"
        }
    },
    PROC_PRES_VAR: {
        "nc_var_name": "pressure",
        "dimension": "time",
        "type": "f8",
        "attrs": {
            "axis": "Z",
            "long_name": "Pressure",
            "observation_type": "measured",
            "positive": "down",
            "reference_datum": "sea-surface",
            "standard_name": "sea_water_pressure",
            "platform": "platform",
            "units": "dbar",
            "instrument": "instrument_ctd",
            "valid_max": 2000.0,
            "valid_min": 0.0,
            "accuracy": 0.01,
            "ancillary_variables": " ",
            "resolution": 0.01,
            "precision": 0.01,
            "_FillValue": nan
        }
    },
    PROC_DEPTH_VAR: {
        "nc_var_name": "depth",
        "dimension": "time",
        "type": "f8",
        "attrs": {
            "_FillValue": nan,
            "axis": "Z",
            "long_name": "Depth",
            "instrument": "instrument_ctd",
            "observation_type": "calculated",
            "positive": "down",
            "reference_datum": "sea-surface",
            "standard_name": "depth",
            "platform": "platform",
            "units": "m",
            "valid_max": 2000.0,
            "valid_min": 0.0,
            "accuracy": 0.01,
            "ancillary_variables": " ",
            "resolution": 0.01,
            "precision": 0.01
        }
    },
    PROC_LAT_VAR: {
        "nc_var_name": "lat",
        "dimension": "time",
        "type": "f8",
        "attrs": {
            "_FillValue": nan,
            "ancillary_variables": " ",
            "coordinate_reference_frame": "urn:ogc:crs:EPSG::4326",
            "reference": "WGS84",
            "platform": "platform",
            "units": "degrees_north",
            "axis": "Y",
            "observation_type": "calculated",
            "long_name": "Estimated Latitude",
            "valid_min": -90.0,
            "valid_max": 90.0,
            "standard_name": "latitude"
        }
    },
    PROC_LON_VAR: {
        "nc_var_name": "lon",
        "dimension": "time",
        "type": "f8",
        "attrs": {
            "_FillValue": nan,
            "ancillary_variables": " ",
            "coordinate_reference_frame": "urn:ogc:crs:EPSG::4326",
            "reference": "WGS84",
            "platform": "platform",
            "units": "degrees_east",
            "axis": "X",
            "observation_type": "calculated",
            "long_name": "Estimated Longitude",
            "valid_min": -180.0,
            "valid_max": 180.0,
            "standard_name": "longitude"
        }
    },
    "sci_water_cond": {
        "nc_var_name": "conductivity",
        "dimension": "time",
        "type": "f8",
        "attrs": {
            "_FillValue": nan,
            "accuracy": 0.0003,
            "ancillary_variables": "conductivity_qc",
            "resolution": 1e-05,
            "precision": "N/A",
            "instrument": "instrument_ctd",
            "long_name": "Conductivity",
            "platform": "platform",
            "observation_type": "measured",
            "standard_name": "sea_water_electrical_conductivity",
            "units": "S m-1",
            "valid_max": 10.0,
            "valid_min": 0.0
        }
    },
    "sci_water_temp": {
        "nc_var_name": "temperature",
        "dimension": "time",
        "type": "f8",
        "attrs": {
            "_FillValue": nan,
            "accuracy": 0.002,
            "ancillary_variables": "temperature_qc",
            "instrument": "instrument_ctd",
            "long_name": "Temperature",
            "observation_type": "measured",
            "platform": "platform",
            "precision": "N/A",
            "resolution": 0.001,
            "standard_name": "sea_water_temperature",
            "units": "Celsius",
            "valid_max": 40.0,
            "valid_min": -5.0
        }
    },
    "salinity": {
        "nc_var_name": "salinity",
        "dimension": "time",
        "type": "f8",
        "attrs": {
            "_FillValue": nan,
            "accuracy": 0.005,
            "ancillary_variables": "salinity_qc",
            "instrument": "instrument_ctd",
            "long_name": "Salinity",
            "observation_type": "calculated",
            "platform": "platform",
            "precision": " ",
            "resolution": " ",
            "standard_name": "sea_water_practical_salinity",
            "units": "1",
            "valid_max": 40.0,
            "valid_min": 0.0
        }
    },
    "density": {
        "nc_var_name": "density",
        "dimension": "time",
        "type": "f8",
        "attrs": {
            "_FillValue": nan,
            "accuracy": " ",
            "ancillary_variables": " ",
            "instrument": "instrument_ctd",
            "long_name": "Density",
            "observation_type": "calculated",
            "platform": "platform",
            "precision": " ",
            "resolution": " ",
            "standard_name": "sea_water_density",
            "units": "kg m-3",
            "valid_max": 1040.0,
            "valid_min": 990.0
        }
    }
}