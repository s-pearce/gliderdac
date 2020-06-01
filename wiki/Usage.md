# How to use gliderdac

gliderdac code is meant to be run over an entire Slocum glider deployment. While single files can be processed, it is meant to run over a directory with glider data, and at this time, particularly Slocum glider ASCII data. I.e. slocum binary data \*bd files that have been converted through the Teledyne Webb Research utilities `dba2asc` and `dba_merge` to a flight and science merged ASCII file.

The main call to a processing run with the code is to the top level python file: [__ooidba_to_ngdac_profile_nc.py__](https://github.com/s-pearce/gliderdac/blob/master/ooidba_to_ngdac_profile_nc.py)

The main and typical inputs are:
* A path to a deployment configuration directory with JSON configuration files (described below)
* The path to the ASCII data files to process as either a list of files to process or a wildcard path that indicates multiple files to process 
    * e.g. data/311/5/*.mrg (I call Slocum glider data that has passed through the `dba_merge` utility, merge files as .mrg)
* An optional path to an output directory where the final netCDF files are created, otherwise they are created in the current directory

So as an example, my typical DAC deployment directory (noted as *) would be
- */config : for the configuration directory and configuration JSON files.
- */nc : for output of the profile NetCDF files
- */logs : for writing any run logs that I capture in the terminal
- */figures : for producing figures to check that profile selection was correct
- and the Slocum ASCII data files (.mrg files) are in a different directory
\<rawdata\>/merged

Then my call would be

```python ooidba_to_ngdac_profile_nc.py -o <dplymnt>/nc <dplymnt>/config <rawdata>/merged/*.mrg >> <dplymnt>/logs/output.log```


## Configuration JSON Files
The deployment configuration directory, of whose path is an input, must have 4 files in it, although it may have additional control files.
If these look familiar it is because they are derived from some of John Kerfoot's Glider Utility controls with some extra abilities.
* [__sensor_defs.json__](https://github.com/s-pearce/gliderdac/blob/master/config-templates/sensor_defs.json) : The definitions that provide variable attributes in the netCDF files, some processing information, and the mapping of Slocum glider data variable names (called "sensors" by TWR) to netCDF variable names.
* [__global_attributes.json__](https://github.com/s-pearce/gliderdac/blob/master/config-templates/global_attributes.json) : The static global attributes that are written to the global attributes in the netCDF file.  Some additional attributes that are specific to the deployment being processed are added by the deployment.json file.  This is meant to provide a static file that doesn't change much that can be copied easily and a shorter one that will change for every deployment processed.
* [__instruments.json__](https://github.com/s-pearce/gliderdac/blob/master/config-templates/instruments.json) : The glider's scientific instrument configuration metadata for the deployment being processed (e.g. Instrument calibration date, serial numbers, make, model, etc.)
* [__deployment.json__](https://github.com/s-pearce/gliderdac/blob/master/config-templates/deployment.json) :  The file that mainly describes the deployment being processed and is unique to the deployment.  It contains information such as a deployment description, WMO ID specific to the glider, deployment launch date, glider name, and other such deployment specific metadata.

For examples of these, see the examples/templates in [config-templates](https://github.com/s-pearce/gliderdac/tree/master/config-templates)
 ToDo: pick a specific directory for these, templates, resources, or examples?

After a deployment is run, a 5th JSON file is created in the configuration directory, __status.json__, which holds information of which input files have already been processed, which profile files have been created, and which have been uploaded to the NGDAC, and if a deployment is a realtime deployment, what is the latest profile ID to continue from for the next run if profile ID is set to be a sequential number (the other option and default is to use the timestamp for the profile).  The __status.json__ file is not really meant to be edited for configuration, but can be done if desired.

An extra feature for these JSON configuration files is that the gliderdac code includes a JSON reader that allows for comments in the JSON files (not a typical JSON syntax). Comments must start with `//`, occupy the whole line (no inline comments), and can include any arbitrary preceding whitespace for indent levels. Since the comments begin with `//` in a JSON file, they are recognized as Java comments by syntax highlighters.
Comments are helpful for remembering options that can be applied in the JSON configuration files. Comments will be added soon to the templates to aid in creation for your own deployments.

__coming soon__: 
As of the current implementation, a copy of each of these files must be present in the configuration directory for the deployment.  Future versions intend to have a master copy of the __sensor_defs.json__ and __global_attributes.json__ files that act as defaults and the local configuration versions can just include elements that need to be changed for the specific deployment (e.g. different comments for a variable).


## Other configuration settings
To aid in control of the running process, there are a few python files that are meant to be changed to meet configuration needs.
* [__configuration.py__](https://github.com/s-pearce/gliderdac/blob/master/configuration.py) : In the top level of the code directory, this is where a list of the science variables from the Slocum glider data should be listed so that the code can check that the relevant science data is present to make sure to include a profile. There are a few other settings that can be made in this file such as some profile filter settings like how many meters deep does a glider need to go, how many minutes in a profile, or how many data points, etc. for a profile to be considered good to keep.
    * In an upcoming version, this will be moved to a different type of configuration file and will be included in a master level configuration directory.
* [__profile_filters.py__](https://github.com/s-pearce/gliderdac/blob/master/profile_filters.py) : This file is meant to be changeable to either change defaults to profile filtering or write your own. To add your own, follow the examples of the other functions and review how a glider data object from this code works [`ooidac.data_classes.GliderData`](https://github.com/s-pearce/gliderdac/blob/master/ooidac/data_classes.py#L15) and prefix the new function you add with `filter_` which the code looks for to run all functions with that prefix.
* [__dba_file_sorter.py__](https://github.com/s-pearce/gliderdac/blob/master/dba_file_sorter.py) : This has a function meant to be changed for sorting the input ASCII data files so that they are run in chronological order.  Chronlogical order is desired so that the Depth Averaged Velocities can be found from the next chronological data file since the correct value for a segment may span several data files.
    * the example included with the code works for OOI filenames that mimic the Slocum binary filenames except with dashes (-) replaced with underscores (_). E.g. ce_311_2019_227_0_19.dbd.asc where ce_311 is the glider name, 2019 the year, 227 the ordinal year day for the start of the mission, 0 the mission number started on ordinal day 227, and 19 being the segment number from mission 227_0.

Glider data processing is highly variable on different operators styles of processing, so of course being an open source code set, the code may be forked and any part edited for your own different use.

The main calling file [__ooidba_to_ngdac_profile_nc.py__](https://github.com/s-pearce/gliderdac/blob/master/ooidba_to_ngdac_profile_nc.py) is the best place to start for editing to suit your own needs, followed by the [processing directory](https://github.com/s-pearce/gliderdac/tree/master/ooidac/processing_dir) and [module](https://github.com/s-pearce/gliderdac/blob/master/ooidac/processing.py), the [profile module](https://github.com/s-pearce/gliderdac/blob/master/ooidac/profiles.py), and the [netCDF writer module](https://github.com/s-pearce/gliderdac/blob/master/ooidac/writers/netCDFwriter.py) under the writers directory.


