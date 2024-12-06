#!/usr/bin/env python
# coding: utf-8
'''
======================================================
Author:  Ömer Özak, 2014 (ozak at smu.edu)
Website: http://omerozak.com
GitHub:  https://github.com/ozak/
======================================================
Package to compute distances using HMI, HMISea and HMIOcean
'''
from __future__ import division
import sys, os, time
from osgeo import gdal, gdalnumeric, ogr, osr
from osgeo.gdalconst import *
import numpy as np
import pandas as pd
from rasterstats import zonal_stats
import pysal as ps
import shapely
from shapely.geometry import Polygon, Point, LineString
from shapely.wkt import loads, dumps
from shapely.ops import cascaded_union
import rasterio.features
from affine import Affine
from shapely.geometry import shape
import georasters as gr
import pyproj
from pyproj import CRS
import geopandas as gp
from geopy.distance import great_circle
import skimage.graph as graph
import multiprocessing as mp
from multiprocessing import Pool, Lock
import warnings
from fiona.crs import from_string
from skimage.segmentation import watershed

path=os.path.dirname(os.path.abspath(__file__))+'/data/'
cea2 = from_string('+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
wgs84 = CRS("EPSG:4326")
cea = CRS("ESRI:54034")

class HMI(object):
    """
    Implements the HMI Class, which will be used to construct the computations, etc. and hold everything in a unique object.

    Usage:
        Load data
        A = HMI(sources, destinations, lat=None, lon=None, x=None, y=None, projected=None, fast=True)
        B = HMI(sources, destinations, lat='lat', lon='lon')
        C = HMI(sources, destinations, y='y', x='x')

        A.correct_all()
        A.HMIdistance(export_shape=True, export_raster=True)
        B.HMIdistance()
        B.hmidist

    where
        sources: (Geo)Pandas DataFrame with the locations of the sources from which to compute distances
                should have columns with either X/Y or LAT/LON
        destinations: (Geo)Pandas DataFrame with the locations of the destinations to which to compute distances
                should have columns with either X/Y or LAT/LON
        lat: Name of column with latitudes (if applicable)
        lon: Name of column with longitudes (if applicable)
        y: Name of column with cylindrical equal area projected latitudes (if applicable)
        x: Name of column with cylindrical equal area projected longitudes (if applicable)
        projected: Boolean, references is data is projected or not
        fast: Use small cost surface or not
        file: Type of cost
        clip: Polygon or GeoPandas DataFrame to clip region for faster local analysis
    """
    # Initialize Object
    def __init__(self, sources, destinations, lat=None, lon=None, x=None, y=None, projected=None, fast=True, file='HMI', clip=None):
        super(HMI, self).__init__()
        self.fast = fast
        self.clip = clip
        self.projected = projected
        self.costs = None
        self.mcp_cost = None
        self.start_pointscols = sources.columns.values
        self.end_pointscols = destinations.columns.values
        if fast:
            self.file = file+'10'
        else:
            self.file = file
        if not isinstance(sources, pd.core.frame.DataFrame) and not isinstance(sources, gp.geodataframe.GeoDataFrame):
            raise TypeError('Sources has to be a (Geo)Pandas Data Frame Object.')
        self.sources = sources.copy()
        if not isinstance(destinations, pd.core.frame.DataFrame) and not isinstance(destinations, gp.geodataframe.GeoDataFrame):
            raise TypeError('Destinations has to be a (Geo)Pandas Data Frame Object.')
        self.destinations = destinations.copy()
        if isinstance(self.destinations, gp.geodataframe.GeoDataFrame):
            projd = self.destinations.crs
        else:
            projd = None
        if isinstance(self.sources, gp.geodataframe.GeoDataFrame):
            projs=self.sources.crs
        else:
            projs = None
        if projd==projs and projd==cea:
            self.projected = True
        else:
            self.projected = False
        if lat!=None and lon!=None:
            self.projected = False
        elif x!=None and y!=None:
            self.projected = True
        if not self.projected:
            warnings.warn('No projection assigned. Program will check Sources and Destinations for clues.\
                          If lat/lon in these data frames, it assumes unprojected. If X/Y in these data frames, then it assumes projected (CEA)')
            lat = self.sources.T.index[self.sources.T.reset_index()['index'].apply(lambda x: x.upper().strip().find('LAT')!=-1)].values
            lon = self.sources.T.index[self.sources.T.reset_index()['index'].apply(lambda x: x.upper().strip().find('LON')!=-1)].values
            latlon = lat.shape[0]+lon.shape[0]
            if latlon>2:
                warnings.warn('More than one columns with words "LAT" or "LON". Chose first occurrence for LAT/LON.\
                                If this is not what you intended, consider assigning `lat` and `lon` strings.')
            if latlon>1:
                self.projected = False
                lat = lat[0]
                lon = lon[0]
            else:
                y = self.sources.T.index[self.sources.T.reset_index()['index'].apply(lambda x: x.upper().strip().find('X')!=-1)].values
                x = self.sources.T.index[self.sources.T.reset_index()['index'].apply(lambda x: x.upper().strip().find('Y')!=-1)].values
                xy = x.shape[0]+y.shape[0]
                if xy>2:
                    warnings.warn('More than one columns with words "X" or "Y". Chose first occurrence for X/Y.\
                                    If this is not what you intended, consider assigning `x` and `y` strings.')
                if xy>1:
                    self.projected =True
                else:
                    raise SyntaxError('Could not determine projection')
        if self.projected==False:
            # Convert project to equal cylindrical area
            p = pyproj.Proj(proj='cea', ellps='WGS84')
            # p(lon,lat)
            [x, y] = p(self.sources[lon].values,self.sources[lat].values)
            self.sources['x'] = x
            self.sources['y'] = y
            [x, y] = p(self.destinations[lon].values,self.destinations[lat].values)
            self.destinations['x'] = x
            self.destinations['y'] = y
        else:
            self.sources['x'] = self.sources[x]
            self.sources['y'] = self.sources[y]
            self.destinations['x'] = self.destinations[x]
            self.destinations['y'] = self.destinations[y]
        self.sources[self.file + 'ID'] = self.sources.index.values
        self.destinations[self.file + 'ID'] = self.destinations.index.values + self.sources[self.file + 'ID'].max() + 1
        self.load_data()
    #############################
    # Set up useful functions
    #############################
    # Load cost surface
    def load_data(self):
        """
        Loads Cost surface as GeoRaster into object's ``costs`` attribute.

        Usage:
            A.load_data()
        """
        if self.fast:
            self.mask = gr.from_file(path+'HMI10.tif')
        else:
            self.mask = gr.from_file(path+'HMI.tif')
        self.costs = gr.from_file(path+self.file+'.tif')
        if self.clip is not None:
            self.costs = self.costs.clip(self.clip)[0]
            self.mask = self.mask.clip(self.clip)[0]
        self.mask = self.mask.raster.mask
        original_mask = self.costs.raster.mask
        self.costs.raster.data[original_mask] = -9
        self.costs.nodata_value = -9
        self.sources['row'], self.sources['col'] = self.costs.map_pixel_location(self.sources.x, self.sources.y)
        self.destinations['row'], self.destinations['col'] = self.costs.map_pixel_location(self.destinations.x, self.destinations.y)
    pass
    # Show cost on HMI cost surface of points at (x,y)
    def costs_points(self, x_points, y_points):
        """
        Show costs at points. Useful to see if any points are outside the costs' surface area.
        Returns a masked Numpy array with the costs

        Usage:
            A.costs_points(A.sources.x, A.sources.y)
        """
        if not self.costs:
            self.load_data()
        return self.costs.map_pixel(x_points, y_points)
    # Costs for sources
    def costs_sources(self):
        """
        Compute costs at sources. Useful to see if any sources are outside the costs' surface area.
        Adds columns to sources with costs

        Usage: A.costs_sources()
        """
        self.sources[self.file+'costs'] = self.costs_points(self.sources.x, self.sources.y)
        self.sources.loc[self.sources[self.file+'costs']==self.costs.nodata_value, self.file+'costs']=np.nan
    pass
    # Costs for destinations
    def costs_destinations(self):
        """
        Compute costs at destinations. Useful to see if any destinations are outside the costs' surface area.
        Adds columns to destinations with costs

        Usage: A.costs_destinations()
        """
        self.destinations[self.file+'costs'] = self.costs_points(self.destinations.x, self.destinations.y)
        self.destinations.loc[self.destinations[self.file+'costs']==self.costs.nodata_value, self.file+'costs']=np.nan
    pass
    # Run both functions
    def costs_all(self):
        """
        Compute costs at both sources and destinations

        Usage: A.costs_all()
        """
        self.costs_sources()
        self.costs_destinations()
    pass
    # Correct bad locations
    def correct_sources(self, maxmove=100, keep=True):
        """
        Move location of sources outside of costs surface up to ``maxmove`` cells in order to have location on cost surface.
        Usage:
            A.correct_destinations(maxmove=100, keep=True)
            B.correct_destinations(maxmove=20, keep=False)
            C.correct_destinations()
        where
            maxmove: Maximum number of cells to move around each location
            keep: Boolean (default True) to create new columns ``originalrow`` and ``originalcol`` in DataFrame to keep original values
        """
        if keep:
            self.sources['originalrow'] = self.sources['row']
            self.sources['originalcol'] = self.sources['col']
        if (self.file+'costs' not in self.sources.columns):
            self.costs_sources()
        problem = self.sources[self.sources[self.file+'costs'].isnull()]
        rows, cols = self.costs.shape
        for i in problem.iterrows():
            xmin = max(i[1].col - maxmove,0)
            xmax = min(i[1].col + maxmove,cols)
            ymin = max(i[1].row - maxmove,0)
            ymax = min(i[1].row + maxmove,rows)
            costn = self.costs[ymin:ymax, xmin:xmax]
            if not isinstance(costn.max(),np.ma.core.MaskedConstant):
                costn = pd.DataFrame(costn>=0).stack().reset_index()
                costn.columns = ['row','col','value']
                costn['row'] = costn['row'] + ymin
                costn['col'] = costn['col'] + xmin
                costn['dist'] = (costn.row-i[1].row)**2 + (costn.col-i[1].col)**2
                costn = costn.loc[costn.dist==costn.dist.min()]
                self.sources.loc[i[0],'row'] = costn.row.values[0]
                self.sources.loc[i[0],'col'] = costn.col.values[0]
    pass
    def correct_destinations(self, maxmove=100, keep=True):
        """
        Move location of destinations outside of costs surface up to ``maxmove`` cells in order to have location on cost surface.
        Usage:
            A.correct_destinations(maxmove=100, keep=True)
            B.correct_destinations(maxmove=20, keep=False)
            C.correct_destinations()
        where
            maxmove: Maximum number of cells to move around each location
            keep: Boolean (default True) to create new columns ``originalrow`` and ``originalcol`` in DataFrame to keep original values
        """
        if keep:
            self.destinations['originalrow'] = self.destinations['row']
            self.destinations['originalcol'] = self.destinations['col']
        if (self.file+'costs' not in self.destinations.columns):
            self.costs_destinations()
        problem = self.destinations[self.destinations[self.file+'costs'].isnull()]
        rows, cols = self.costs.shape
        for i in problem.iterrows():
            xmin = max(i[1].col - maxmove,0)
            xmax = min(i[1].col + maxmove,cols)
            ymin = max(i[1].row - maxmove,0)
            ymax = min(i[1].row + maxmove,rows)
            costn = self.costs[ymin:ymax, xmin:xmax]
            if not isinstance(costn.max(),np.ma.core.MaskedConstant):
                costn = pd.DataFrame(costn>=0).stack().reset_index()
                costn.columns = ['row','col','value']
                costn['row'] = costn['row'] + ymin
                costn['col'] = costn['col'] + xmin
                costn['dist'] = (costn.row-i[1].row)**2 + (costn.col-i[1].col)**2
                costn = costn.loc[costn.dist==costn.dist.min()]
                self.destinations.loc[i[0],'row'] = costn.row.values[0]
                self.destinations.loc[i[0],'col'] = costn.col.values[0]
    pass
    def correct_data(self, maxmove=100, keep=True):
        """
        Move location of sources and destinations outside of costs surface up to ``maxmove`` cells in order to have location on cost surface.
        Usage:
            A.correct_data(maxmove=100, keep=True)
            B.correct_data(maxmove=20, keep=False)
            C.correct_data()
        where
            maxmove: Maximum number of cells to move around each location
            keep: Boolean (default True) to create new columns ``originalrow`` and ``originalcol`` in DataFrame to keep original values
        """
        self.correct_sources(maxmove=maxmove, keep=keep)
        self.correct_destinations(maxmove=maxmove, keep=keep)
    pass
    # Setup Graph for distance computations and provide distance functions
    def mcp(self, *args, **kwargs):
        """
        Setup MCP_Geometric object from skimage for optimal travel time computations
        """
        # Create Cost surface to work on
        self.mcp_cost=graph.MCP_Geometric(self.costs.raster, *args, **kwargs)
    pass
    # Determine minimum travel cost to each location
    def HMIdistance(self, isolation=True, export_raster=False, export_shape=False, routes=False, path='./'):
        """
        Compute HMI based travel distance measured in weeks of travel from each start point to all end points.
        The function returns the distances between the start point and the end points as a Pandas dataframe.
        Additionally, for each start point it computes the level of isolation, i.e. its average travel distance to
        all other location on land
        """
        start_points=self.sources
        end_points=self.destinations
        if not self.mcp_cost:
            self.mcp()
        count=0
        for i in start_points.iterrows():
            cumulative_costs, traceback = self.mcp_cost.find_costs([[i[1].row,i[1].col]])
            dist=cumulative_costs[end_points.row.values,end_points.col.values].transpose()/(7*24)
            df2=pd.DataFrame(np.array([(i[1][self.file+'ID']*np.ones_like(dist)).flatten(),
                                    end_points[self.file+'ID'],dist.flatten()]).transpose(),
                                    columns=[self.file+'ID1',self.file+'ID2',self.file+'dist'])
            # Keep only locations that are accessible
            df2 = df2.loc[df2[self.file+'dist']<np.inf]
            if isolation:
                hmiisolation = np.ma.masked_array(cumulative_costs, mask = np.logical_or(self.mask, cumulative_costs==np.inf)
                                    , fill_value=np.nan).mean()/(7*24)
                start_points.loc[i[0],self.file+'Iso'] = hmiisolation
            if export_raster:
                cumulative_costs = gr.GeoRaster(np.ma.masked_array(cumulative_costs,
                                                    mask = np.logical_or(self.costs.raster.mask, cumulative_costs==np.inf),
                                                    fill_value = np.nan), self.costs.geot, self.costs.nodata_value,
                                                    projection=self.costs.projection, datatype = self.costs.datatype)
                cumulative_costs.raster.data[cumulative_costs.raster.mask] = cumulative_costs.nodata_value
                cumulative_costs.to_tiff(path+self.file+'_ID'+str(i[1][self.file+'ID']))
            if df2.size>0:
                if export_shape:
                    routes=True
                if routes:
                    df2['geometry'] = df2[self.file+'ID2'].apply(lambda x:
                                                self.mcp_cost.traceback(end_points.loc[end_points[self.file+'ID']==x][['row','col']].values[0]))
                    df2['geometry'] = df2.geometry.apply(lambda x: [gr.map_pixel_inv(y[0],y[1], self.costs.geot[1],
                                                    self.costs.geot[-1], self.costs.geot[0], self.costs.geot[-3]) for y in x])
                    df2['geometry'] = df2.geometry.apply(lambda x: LineString(x) if int(len(x)>1) else LineString([x[0],x[0]]))
                    df2 = gp.GeoDataFrame(df2, crs = cea)
                if isolation:
                    df2[self.file+'Iso'] = hmiisolation
                if count==0:
                    self.hmidist = df2.copy()
                else:
                    self.hmidist = pd.concat([self.hmidist, df2])
                count+=1
        if routes:
            self.hmidist = gp.GeoDataFrame(self.hmidist, crs = cea)
        if export_shape:
            if 'geometry' in self.end_pointscols:
                self.hmidist = pd.merge(self.hmidist, end_points[[self.file+'ID']+self.end_pointscols.tolist()].drop('geometry', axis=1), left_on=self.file+'ID2', right_on=self.file+'ID', how='left')
            else:
                self.hmidist = pd.merge(self.hmidist, end_points[[self.file+'ID']+self.end_pointscols.tolist()], left_on=self.file+'ID2', right_on=self.file+'ID', how='left')
            if 'geometry' in self.start_pointscols:
                self.hmidist = pd.merge(self.hmidist, start_points[[self.file+'ID']+self.start_pointscols.tolist()].drop('geometry', axis=1), left_on=self.file+'ID1', right_on=self.file+'ID', how='left',
                             suffixes = ['_2','_1'])
            else:
                self.hmidist = pd.merge(self.hmidist, start_points[[self.file+'ID']+self.start_pointscols.tolist()], left_on=self.file+'ID1', right_on=self.file+'ID', how='left',
                             suffixes = ['_2','_1'])
            self.hmidist = gp.GeoDataFrame(self.hmidist, crs = cea)
            self.hmidist.to_file(path+self.file+'_routes.geojson', driver="GeoJSON")
    pass
    # Create Voronoi
    def voronoi(self):
        if not self.mcp_cost:
            self.mcp()
        mat_sources = np.ones(self.costs.shape, dtype=int) * 0
        for i in self.sources.iterrows():
            mat_sources[i[1].row, i[1].col] = i[1][self.file + 'ID'] + 1
        cumulative_costs, traceback = self.mcp_cost.find_costs(self.sources[['row', 'col']].values)
        mask = np.less(cumulative_costs, np.inf)
        allo = watershed(cumulative_costs, mat_sources, mask=mask)
        allo = allo - 1
        allo[allo==-1] = self.costs.nodata_value
        voronoi_raster = gr.GeoRaster(allo, self.costs.geot, nodata_value=self.costs.nodata_value, projection=self.costs.projection, datatype=self.costs.datatype)
        voronoi = []
        for shp, val in rasterio.features.shapes(allo, transform=Affine.from_gdal(*self.costs.geot)):
            voronoi.append([val, shape(shp)])
        voronoi = gp.GeoDataFrame(voronoi, crs=cea, columns=[self.file + 'ID', 'geometry'])
        voronoi = voronoi.merge(self.sources[self.sources.columns.difference(['geometry'])], on=self.file + 'ID')
        return voronoi

# Implement HMISea class, whcih is basically the same except for the raster being used
class HMISea(HMI):
    """
    Implements the HMISea Class, which will be used to construct the computations, etc. and hold everything in a unique object.

    Usage:
        Load data
        A = HMISea(sources, destinations, lat=None, lon=None, x=None, y=None, projected=None, fast=True)
        B = HMISea(sources, destinations, lat='lat', lon='lon')
        C = HMISea(sources, destinations, y='y', x='x')

        A.correct_all()
        A.HMIdistance(export_shape=True, export_raster=True)
        B.HMIdistance()
        B.hmidist

    where
        sources: (Geo)Pandas DataFrame with the locations of the sources from which to compute distances
                should have columns with either X/Y or LAT/LON
        destinations: (Geo)Pandas DataFrame with the locations of the destinations to which to compute distances
                should have columns with either X/Y or LAT/LON
        lat: Name of column with latitudes (if applicable)
        lon: Name of column with longitudes (if applicable)
        y: Name of column with cylindrical equal area projected latitudes (if applicable)
        x: Name of column with cylindrical equal area projected longitudes (if applicable)
        projected: Boolean, references is data is projected or not
        fast: Use small cost surface or not
    """
    # Initialize Object
    def __init__(self, sources, destinations, lat=None, lon=None, x=None, y=None, projected=None, fast=True, file='HMISea', clip=None):
        super(HMISea, self).__init__(sources, destinations, lat=lat, lon=lon, x=x, y=y, projected=projected, fast=fast, file='HMISea', clip=clip)
    pass
    def HMISeadistance(self, isolation=True, export_raster=False, export_shape=False, routes=False, path='./'):
        """
        Compute HMI based travel distance measured in weeks of travel from each start point to all end points.
        The function returns the distances between the start point and the end points as a Pandas dataframe.
        Additionally, for each start point it computes the level of isolation, i.e. its average travel distance to
        all other location on land
        """
        self.HMIdistance(isolation=isolation, export_raster=export_raster, export_shape=export_shape, routes=routes, path=path)
        self.hmiseadist = self.hmidist
    pass

'''
# Implement HMIOcean class, whcih is basically the same except for the raster being used
class HMIOcean(HMI):
    """
    Implements the HMISea Class, which will be used to construct the computations, etc. and hold everything in a unique object.

    Usage:
        Load data
        A = HMISea(sources, destinations, lat=None, lon=None, x=None, y=None, projected=None, fast=True)
        B = HMISea(sources, destinations, lat='lat', lon='lon')
        C = HMISea(sources, destinations, y='y', x='x')

        A.correct_all()
        A.HMIdistance(export_shape=True, export_raster=True)
        B.HMIdistance()
        B.hmidist

    where
        sources: (Geo)Pandas DataFrame with the locations of the sources from which to compute distances
                should have columns with either X/Y or LAT/LON
        destinations: (Geo)Pandas DataFrame with the locations of the destinations to which to compute distances
                should have columns with either X/Y or LAT/LON
        lat: Name of column with latitudes (if applicable)
        lon: Name of column with longitudes (if applicable)
        y: Name of column with cylindrical equal area projected latitudes (if applicable)
        x: Name of column with cylindrical equal area projected longitudes (if applicable)
        projected: Boolean, references is data is projected or not
        fast: Use small cost surface or not
    """
    # Initialize Object
    def __init__(self, sources, destinations, lat=None, lon=None, x=None, y=None, projected=None, fast=True, file='HMIOcean', clip=None):
        super(HMIOcean, self).__init__(sources, destinations, lat=lat, lon=lon, x=x, y=y, projected=projected, fast=fast, file='HMIOcean', clip=clip)
    pass
    def HMIOceandistance(self, isolation=True, export_raster=False, export_shape=False, routes=False, path='./'):
        """
        Compute HMI based travel distance measured in weeks of travel from each start point to all end points.
        The function returns the distances between the start point and the end points as a Pandas dataframe.
        Additionally, for each start point it computes the level of isolation, i.e. its average travel distance to
        all other location on land
        """
        self.HMIdistance(isolation=isolation, export_raster=export_raster, export_shape=export_shape, routes=routes, path=path)
        self.hmioceandist = self.hmidist
    pass
'''

def download():
    path = os.path.dirname(os.path.realpath(__file__)) + '/data/'
    urlretrieve('https://zenodo.org/records/14285746/files/HMI.tif?download=1', filename=path+'HMI.tif')
    urlretrieve('https://zenodo.org/records/14285746/files/HMISea.tif?download=1', filename=path+'HMISea.tif')
    urlretrieve('https://zenodo.org/records/14285746/files/HMI10.tif?download=1', filename=path+'HMI10.tif')
    urlretrieve('https://zenodo.org/records/14285746/files/HMISea10.tif?download=1', filename=path+'HMISea10.tif')    
    pass
