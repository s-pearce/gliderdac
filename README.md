# gliderdac - IOOS Glider DAC file creation
### Warning: This repository is still under development!

The U.S. [Integrated Ocean Observing System \(IOOS\)](https://gliders.ioos.us/) hosts the [National Glider Data Assembly \(DAC\)](https://gliders.ioos.us/data) for scientific data collected by gliders (autonomous underwater vehicles).  `gliderdac` is a python package for writing IOOS Glider DAC complient NetCDF files from [Ocean Observatories Initiative \(OOI\)](https://oceanobservatories.org/) specific [Slocum](http://www.teledynemarine.com/slocum-glider) glider data files.  `gliderdac` was adapted from John Kerfoot's
Glider NetCDF Utilities [gncutils](https://github.com/kerfoot/gncutils)

Full documentation available at the [wiki](https://github.com/s-pearce/gliderdac/wiki)

## Requirements
+ Python 3.*
+ [netCDF4](https://unidata.github.io/netcdf4-python/netCDF4/index.html)
+ [numpy](https://www.numpy.org/)
+ [scipy](https://www.scipy.org/)
+ [matplotlib](https://matplotlib.org/)
+ [GSW-Python](https://github.com/TEOS-10/GSW-Python)
+ [Shapely](https://pypi.org/project/Shapely/)

## Installation

    > git clone https://github.com/s-pearce/gliderdac.git
    > cd gliderdac
    > conda env create -f environment.yml



## License
GNU General Public License v3.0

See LICENSE.txt to see the full text.
