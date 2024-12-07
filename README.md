The Human Mobility Index Python package - `hmi`
===========

<a href="https://pypi.python.org/pypi/hmi/">![PyPiVersion](https://img.shields.io/pypi/v/hmi.svg)</a> <a href="">![Pyversions](https://img.shields.io/pypi/pyversions/hmi.svg)</a> <a href="https://hmi.readthedocs.io/en/latest/">![ReadTheDocs](https://readthedocs.org/projects/hmi/badge/?version=latest&style=plastic)</a> <a>[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14291903.svg)](https://doi.org/10.5281/zenodo.14291903)</a>

The ``hmi`` package is a python module that provides an interface to compute distances between locations based on [Özak's (2018, ](http://rdcu.be/I4YI)[2010)](http://omerozak.com/pdf/Ozak_voyage.pdf) [Human Mobility Index](https://human-mobility-index.github.io/). There are three main classes of computations based on the Human Mobility Index:

- HMI: Only land based mobility
- HMISea: Land based mobility and pre-1500CE seafaring technologies
- HMIOcean: Land based mobility and post-1500CE seafaring technologies

It includes tools to 

- Compute distances between locations
- Compute isolation measures for locations (see e.g. Ashraf, Galor, and Özak (2010))
- Generate optimal routes as shape files
- Export optimal travel times, routes, etc. to shape files, CSV files, etc

Install
-------

Given the requirements of the package, the best way to install is to first create a `mamba/conda`  environment as follows (this creates an `python-3.11` environment and adds some basic packages that are needed).

```bash
mamba create --name HMI --override-channels -c conda-forge python=3.11 pip geopandas georasters jupyterlab jupyter seaborn geoplot pysal 
```

activate the environment

```bash
mamba activate HMI
```

then 

```bash
 pip install hmi
```
 or

```bash
 pip install git+git://github.com/ozak/hmi.git
```

Example Usage: HMI
-------------------------

``` python
 import hmi 
``` 

## Citation

If you use the package please cite:

[Özak, Ömer. "Distance to the pre-industrial technological frontier and economic development." Journal of Economic Growth 23.2 (2018): 175-221.](http://rdcu.be/I4YI)

[Özak, Ömer. "The voyage of homo-economicus: Some economic measures of distance." Department of Economics, Brown University (2010).](http://omerozak.com/pdf/Ozak_voyage.pdf)

and the package itself.
 
Issues
------

Find a bug? Report it via github issues by providing

- a link to download the smallest possible raster and vector dataset necessary to reproduce the error
- python code or command to reproduce the error
- information on your environment: versions of python, gdal and numpy and system memory


# Copyright 

&copy; Ömer Özak (2014)

This code and data is provided under [Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) License](https://creativecommons.org/licenses/by-sa/4.0/) and [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html).
![](http://mirrors.creativecommons.org/presskit/buttons/88x31/svg/by-sa.svg) ![](https://www.gnu.org/graphics/gplv3-127x51.png)
