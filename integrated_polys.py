#!/usr/bin/env python

__author__ = "Natasha Hurley-Walker"
__date__ = "11/09/2018"

import os
import sys
import shutil
import glob

import matplotlib
matplotlib.use('Agg') # So does not use display -- only good if just making plots
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse

from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['serif']})

from astropy.io import fits
from astropy import wcs
#from astropy.table import Table, Column
from astropy import units as u
from astropy.nddata import Cutout2D

from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve

import numpy as np

import defs
from defs import SNR
from defs import read_snrs
from defs import find_fluxes
from measure_snr import import_snr
from measure_source import make_single_plot

def Jy2mJy(S):
    return S*1000

def mJy2Jy(S):
    return S/1000

def deg2AS(angle):
    return angle*3600

def WL2freq(WL):
    return 3.e8/WL

def freq2WL(nu):
    return 3.e8/nu

def m2cm(length):
    return length*100

def cm2m(length):
    return length/100

def mK2K(T):
    return T/1000

#https://science.nrao.edu/facilities/vla/proposing/TBconv
def T2PFD(T, nu, bmaj, bmin):
    ''' Convert a brightness temperature to peak flux density. \
        Expecting T in K, nu in Hz, bmaj and bmin in deg; \
        Will convert to mJy, cm, arcsec." '''
    PFD = T * deg2AS(bmaj) * deg2AS(bmin) /     \
       ( 1.36 * (m2cm(freq2WL(nu)))**2 )
    return mJy2Jy(PFD)

def measure_E11(snr):
    sourcedir = "/media/data/Dropbox/MWA/SNR_search/confirm_candidates/regrid_all/E11/"
    fitsfile = sourcedir+snr.name+"_E11.fits"
    final_flux, total_flux, bkg_flux, rms_flux, nbeam = find_fluxes(snr.polygon, snr.sources, snr.exclude, fitsfile)
    hdu = fits.open(fitsfile)
    T = mK2K(final_flux)
    nu = WL2freq(cm2m(11.)) # Effelsberg 11cm survey
    bmaj = hdu[0].header["BMAJ"]
    bmin = hdu[0].header["BMIN"]
    make_single_plot(snr.polygon, snr.sources, snr.exclude, fitsfile)
    print "final flux, total flux, background flux, rms flux, number of beams"
    print final_flux, total_flux, bkg_flux, rms_flux, nbeam
    S_11 = T2PFD(T, nu, bmaj, bmin)
    return S_11

if __name__ == "__main__":
    snrs = read_snrs()
    updated_snrs = []
    for snr in snrs:
        snr = import_snr(snr)
        updated_snrs.append(snr)
    snrs = updated_snrs

    for snr in snrs:
        S_11 = measure_E11(snr)
        print S_11
