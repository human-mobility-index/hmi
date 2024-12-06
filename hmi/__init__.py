# -*- coding: utf-8 -*-
"""The HMI package is a python module that provides an interface to compute distances between locations based on Özak's (2010, 2012, 2015) Human Mobility Index."""

#import metadata
from .metadata import *
from .hmi import wgs84, cea, HMI, HMISea, path
import os
import requests
from tqdm import tqdm

__version__ = version
__author__ = authors[0]
__license__ = license
__copyright__ = copyright

def download_file_with_progress(url, filename):
    """Download a file with a progress bar."""
    response = requests.get(url, stream=True)
    total_size = int(response.headers.get('content-length', 0))
    block_size = 1024  # 1 Kibibyte

    # Use tqdm for the progress bar
    with tqdm(total=total_size, unit='iB', unit_scale=True, desc=f"Downloading {filename}") as t:
        with open(filename, 'wb') as f:
            for data in response.iter_content(block_size):
                t.update(len(data))
                f.write(data)

# Download data if needed
def download():
    path = os.path.dirname(os.path.realpath(__file__)) + '/data/'
    os.makedirs(path, exist_ok=True)

    files = {
        'HMI.tif': 'https://zenodo.org/records/14285746/files/HMI.tif?download=1',
        'HMISea.tif': 'https://zenodo.org/records/14285746/files/HMISea.tif?download=1',
        'HMI10.tif': 'https://zenodo.org/records/14285746/files/HMI10.tif?download=1',
        'HMISea10.tif': 'https://zenodo.org/records/14285746/files/HMISea10.tif?download=1',
    }

    for filename, url in files.items():
        dest_path = os.path.join(path, filename)
        if not os.path.exists(dest_path):
            download_file_with_progress(url, dest_path)

# Check whether data is present or not
path = os.path.dirname(os.path.realpath(__file__)) + '/data/'
if not os.path.exists(path + 'HMI10.tif'):
    print('Downloading required HMI datasets.')
    print('This may take a few minutes.')
    print('This will be done only once.')
    print('Please be patient.')
    download()