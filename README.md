# gliderdac - IOOS Glider DAC file creation
# This repo is still under development

The U.S. Integrated Ocean Observing System \([IOOS](https://gliders.ioos.us/)\) hosts the [National Glider Data Assembly] (https://gliders.ioos.us/data/) \(DAC\) for scientific data collected by gliders (autonomous underwater vehicles).  gliderdac is a python package for writing IOOS [GliderDAC] complient NetCDF files from [The Ocean Observatories Initiative](https://oceanobservatories.org/) \(OOI\) specific Slocum glider data files.  Loosely based off of John Kerfoot's
Glider NetCDF Utilities [wiki](https://github.com/kerfoot/gncutils/wiki)

Full documentation available at the [wiki](https://github.com/s-pearce/gliderdac/wiki)

# Requirements
+ Python 3.*
+ Conda
+ [netCDF4](https://unidata.github.io/netcdf4-python/netCDF4/index.html)
+ [numpy](https://www.numpy.org/)
+ [scipy](https://www.scipy.org/)

# Installation

    > git clone https://github.com/s-pearce/gliderdac.git
    > cd gliderdac
    > conda env create -f environment.yml



# License
GNU General Public License v3.0

See LICENSE.txt to see the full text.
