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
