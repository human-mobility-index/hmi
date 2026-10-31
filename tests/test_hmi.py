# coding: utf-8
'''
End-to-end tests for the hmi package using the bundled cost surfaces (fast=True).
Run from the repository root with: pytest
'''
import warnings
import numpy as np
import pandas as pd
import pytest

warnings.simplefilter('ignore')
import hmi

# Sources get IDs 0-1, destinations 2-5 (in this order)
SOURCES = pd.DataFrame({'name': ['Addis Ababa', 'London'],
                        'lat': [9.03, 51.51], 'lon': [38.74, -0.13]})
DESTINATIONS = pd.DataFrame({'name': ['Cairo', 'Beijing', 'Lima', 'Sydney'],
                             'lat': [30.04, 39.90, -12.05, -33.87],
                             'lon': [31.24, 116.41, -77.04, 151.21]})


def run(cls):
    A = cls(SOURCES, DESTINATIONS, lat='lat', lon='lon')
    A.costs_all()
    A.correct_data()
    A.HMIdistance(routes=True)
    d = A.hmidist.set_index([A.file + 'ID1', A.file + 'ID2'])
    return A, d


@pytest.fixture(scope='module')
def hmi_land():
    return run(hmi.HMI)


@pytest.fixture(scope='module')
def hmi_sea():
    return run(hmi.HMISea)


def test_download_lives_in_init():
    assert callable(hmi.download)
    assert not hasattr(hmi.hmi, 'download')


def test_projection(hmi_land):
    A, _ = hmi_land
    assert not A.projected
    assert {'x', 'y', 'row', 'col'} <= set(A.sources.columns)


def test_hmi_land_distances(hmi_land):
    A, d = hmi_land
    dist = d[A.file + 'dist']
    # Only Addis Ababa reaches Cairo and Beijing over land
    assert set(d.index) == {(0, 2), (0, 3)}
    assert dist[(0, 2)] == pytest.approx(4.08, abs=0.01)
    assert dist[(0, 3)] == pytest.approx(12.36, abs=0.01)
    assert np.isfinite(d[A.file + 'Iso']).all()
    assert d.geometry.notna().all()


def test_hmi_sea_distances(hmi_sea, hmi_land):
    A, d = hmi_sea
    dist = d[A.file + 'dist']
    # Sea travel adds London and Sydney, but not Lima (Americas)
    assert set(d.index) == {(s, t) for s in (0, 1) for t in (2, 3, 5)}
    assert dist[(0, 2)] == pytest.approx(3.32, abs=0.01)
    assert dist[(1, 5)] == pytest.approx(19.80, abs=0.01)
    # Sea access can only shorten trips
    Al, dl = hmi_land
    assert dist[(0, 3)] < dl[Al.file + 'dist'][(0, 3)]


def test_voronoi():
    v = hmi.HMI(SOURCES, DESTINATIONS, lat='lat', lon='lon').voronoi()
    assert v.crs == hmi.cea
    assert sorted(v.name.unique()) == ['Addis Ababa', 'London']
    assert v.geometry.is_valid.all()


# Points just off the coast (in the sea) that correct_data should move onto land
OFFSHORE = pd.DataFrame({'name': ['Cairo', 'Off Alexandria', 'Off Beirut'],
                         'lat': [30.04, 31.4, 33.9], 'lon': [31.24, 29.9, 35.3]})


def test_correct_data_moves_offshore_points():
    A = hmi.HMI(SOURCES.iloc[:1], OFFSHORE, lat='lat', lon='lon')
    A.costs_all()
    off = A.destinations[A.file + 'costs'].isnull()
    assert off.tolist() == [False, True, True]
    A.correct_data(maxmove=10)
    moved = A.destinations.loc[off]
    # Moved onto valid cells, within maxmove cells (Euclidean) of the original location
    assert (A.costs.raster[moved.row.values, moved.col.values] >= 0).all()
    assert not np.ma.is_masked(A.costs.raster[moved.row.values, moved.col.values])
    d2 = (moved.row - moved.originalrow)**2 + (moved.col - moved.originalcol)**2
    assert ((d2 > 0) & (d2 <= 10**2)).all()
    # Points on the surface stay put
    assert (A.destinations.loc[~off, 'row'] == A.destinations.loc[~off, 'originalrow']).all()
    # All destinations, including the last one, are now reachable
    A.HMIdistance()
    assert sorted(A.hmidist[A.file + 'ID2'].astype(int)) == [1, 2, 3]


def test_correct_data_respects_maxmove():
    A = hmi.HMI(SOURCES.iloc[:1], OFFSHORE, lat='lat', lon='lon')
    A.correct_data(maxmove=0)
    assert (A.destinations.row == A.destinations.originalrow).all()
    assert (A.destinations.col == A.destinations.originalcol).all()
    A.HMIdistance()
    assert sorted(A.hmidist[A.file + 'ID2'].astype(int)) == [1]


def test_floor_option():
    A = hmi.HMI(SOURCES, DESTINATIONS, lat='lat', lon='lon')
    B = hmi.HMI(SOURCES, DESTINATIONS, lat='lat', lon='lon', floor=True)
    assert not A.floor and B.floor
    g = A.costs
    for df in (A.destinations, A.sources):
        col = (df.x - g.xmin) / g.x_cell_size
        row = (df.y - g.ymax) / g.y_cell_size
        assert (df.col == np.round(col)).all() and (df.row == np.round(row)).all()
    for df in (B.destinations, B.sources):
        col = (df.x - g.xmin) / g.x_cell_size
        row = (df.y - g.ymax) / g.y_cell_size
        assert (df.col == np.floor(col)).all() and (df.row == np.floor(row)).all()
    # costs_points uses the same rule as the stored rows/cols
    B.costs_destinations()
    expected = B.costs.raster[B.destinations.row.values, B.destinations.col.values]
    assert np.allclose(B.destinations[B.file + 'costs'], expected)
    assert hmi.HMISea(SOURCES, DESTINATIONS, lat='lat', lon='lon', floor=True).floor


def test_correct_ties_are_deterministic():
    import georasters as gr
    A = hmi.HMI(SOURCES.iloc[:1], DESTINATIONS.iloc[:1], lat='lat', lon='lon')
    # 5x5 surface: only the four cells at distance 1 from the centre (2,2) are valid
    data = np.full((5, 5), -9.0)
    for r, c in [(1, 2), (2, 1), (2, 3), (3, 2)]:
        data[r, c] = 1.0
    A.costs = gr.GeoRaster(np.ma.masked_equal(data, -9.0), (0, 1, 0, 0, 0, -1), nodata_value=-9)
    pts = pd.DataFrame({'row': [2], 'col': [2], A.file + 'costs': [np.nan]})
    A._correct_points(pts, maxmove=2)
    # Tie among the four neighbours goes to the northernmost (smallest row)
    assert (pts.row[0], pts.col[0]) == (1, 2)
    # Without the northern cell, the tie goes to the westernmost of the remaining row
    data[1, 2] = -9.0
    A.costs = gr.GeoRaster(np.ma.masked_equal(data, -9.0), (0, 1, 0, 0, 0, -1), nodata_value=-9)
    pts = pd.DataFrame({'row': [2], 'col': [2], A.file + 'costs': [np.nan]})
    A._correct_points(pts, maxmove=2)
    assert (pts.row[0], pts.col[0]) == (2, 1)
