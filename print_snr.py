#!/usr/bin/env python

__author__ = "Natasha Hurley-Walker"
__date__ = "20/12/2018"

import os
import sys
import shutil
import glob

import numpy as np
import pickle
from astropy.io import fits
from astropy import wcs
from astropy.table import Table, Column

import argparse

import defs
from defs import Coords
from defs import SNR
from defs import read_snrs

snrs = read_snrs()

for snr in snrs:
    print snr.name

