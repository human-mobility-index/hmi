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

The package depends on GDAL (through `georasters` and `rasterio`), which `pip` cannot build without a system GDAL library. The best way to install is to first create a `mamba/conda` environment with the geospatial dependencies from `conda-forge`:

```bash
mamba create --name HMI --override-channels -c conda-forge python=3.12 pip gdal geopandas georasters rasterio scikit-image requests tqdm
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
pip install git+https://github.com/human-mobility-index/hmi.git
```

The first time you `import hmi`, the package downloads the HMI cost surfaces (about 4.5 GB) from [Zenodo](https://doi.org/10.5281/zenodo.14285746) into its `data/` folder. This is done only once.

Example Usage: HMI
-------------------------

Sources and destinations are (Geo)Pandas data frames with latitude and longitude (or projected CEA `x`/`y`) columns. Distances are measured in weeks of travel.

``` python
import pandas as pd
import hmi

sources = pd.DataFrame({'name': ['Addis Ababa', 'London'],
                        'lat': [9.03, 51.51], 'lon': [38.74, -0.13]})
destinations = pd.DataFrame({'name': ['Cairo', 'Beijing', 'Sydney'],
                             'lat': [30.04, 39.90, -33.87], 'lon': [31.24, 116.41, 151.21]})

# Land-only mobility (use hmi.HMISea to add pre-1500CE seafaring)
A = hmi.HMI(sources, destinations, lat='lat', lon='lon')
A.correct_data()             # move points that fall off the cost surface onto it
A.HMIdistance(routes=True)   # travel times, isolation and optimal routes
A.hmidist                    # GeoDataFrame with one row per reachable pair

# Regions closest (in travel time) to each source
regions = A.voronoi()
```

`correct_data(maxmove=100)` moves each point that falls off the cost surface (e.g. a coastal city whose cell is in the sea) to the nearest valid cell within `maxmove` cells; points with no valid cell that close are left where they are. By default points are assigned to cells by rounding their coordinates, as in `georasters`; pass `floor=True` to `hmi.HMI`/`hmi.HMISea` to use the cell that contains each point instead.

Pairs that cannot be reached (e.g. London over land only, or points that are still off the cost surface) are dropped from `hmidist`. Use `A.HMIdistance(export_shape=True, export_raster=True, path='./')` to write the routes as GeoJSON and the travel-time surfaces as GeoTIFFs.

Testing
-------

From a clone of the repository, inside the environment above:

```bash
pip install pytest
pytest
```

The tests compute HMI and HMISea distances, isolation, routes and Voronoi regions for a few world cities and compare them with known values.

## Citation

If you use the package please cite:

[Özak, Ömer. "Distance to the pre-industrial technological frontier and economic development." Journal of Economic Growth 23.2 (2018): 175-221.](http://rdcu.be/I4YI)

[Özak, Ömer. "The voyage of homo-economicus: Some economic measures of distance." Department of Economics, Brown University (2010).](http://omerozak.com/pdf/Ozak_voyage.pdf)

and the package itself

[Özak, Ömer. “The Human Mobility Index Python Package - `hmi`”. Zenodo, December 7, 2024. https://doi.org/10.5281/zenodo.14291912.](https://doi.org/10.5281/zenodo.14291912)
 
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
