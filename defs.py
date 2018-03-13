#!/usr/bin/env python

__author__ = "Natasha Hurley-Walker"
__date__ = "10/03/2018"

import os
import sys
from astropy.io import fits
#from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u

import numpy as np

# SNR naming format
global nameformat
nameformat="SNR_G{0:003.1f}{1:+0.1f}"

# Line markers and colors
firstsnr = { "marker" : "x" , "linestyle" : "None", "color" : "red" }
restsnr = { "marker" : "None" , "linestyle" : "-", "color" : "red" }
srcmark = { "marker" : "o" , "linestyle" : "None", "color" : "yellow" }
firstexc = { "marker" : "x" , "linestyle" : "None", "color" : "blue" }
restexc = { "marker" : "None" , "linestyle" : ":", "color" : "blue" }
reticsnr = { "marker" : "None" , "linestyle" : "--", "color" : "green" }

class Coords:
    def __init__(self):
        self.x = []
        self.y = []

class SNR(object):
    def __init__(self):
        self.name = None
        self.loc = None
        self.maj = None
        self.min = None
        self.known = None
        self.polygon = None
        self.sources = None
        self.exclude = None
        self.flux = None # A dictionary of frequencies to flux densities
        self.fit = None # A dictionary of "flux" -> total flux at 150 MHz; "alpha" -> fitted alpha; "chi2red" -> reduced chi2

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
# Read in the SNRs from files
# Get the known SNRs from Dave Green's "DAG.txt"
#                                          l b HH MM SS DD MM 
def read_snrs():
    if not os.path.exists("DAG.txt") or not os.path.exists("candidates.txt"):
        print "Argh, can't find DAG.txt and candidates.txt!"
        sys.exit(1)
        
    snrs = []

    if file_len("DAG.txt"):
        tempcat = np.genfromtxt("DAG.txt", delimiter=(6,6,4,3,3,5,3))
        # Elliptical radii
        sizestrs = np.genfromtxt("DAG.txt",usecols=7,dtype=str)
# Makes it work for single-line files
        ac = 1
        if not len(sizestrs.shape):
            ac = 0
            sizestrs = np.array([sizestrs])
        a = []
        b = []
        for sizestr in sizestrs:
           a.append(sizestr.replace("?","").split("x")[0])
           if "x" in sizestr:
               b.append(sizestr.replace("?","").split("x")[1])
           else:
                b.append(sizestr.replace("?",""))
        snrcat = np.concatenate([tempcat,np.squeeze(np.array([a,b],dtype='float')).transpose()],axis=ac)
# Makes it work for single-line files
        if len(snrcat.shape)==1:
           snrcat = [snrcat]

        # Transform the array into a list of SNR objects
        for row in snrcat:
           snr = SNR()
        #   l = row[0]
        #   b = row[1]
        #   RA = 15*(row[2]+(row[3]/60.)+(row[4]/3600.))
        #   Dec = row[5]+(row[6]/60.)
           snr.loc = SkyCoord("{0:02.0f}h{1:02.0f}m{2:02.0f}s {3:+02.0f}d{4:02.0f}m00s".format(row[2],row[3],row[4],row[5],row[6]))
           snr.maj = row[7]/60.
           snr.min = row[8]/60.
           snr.name  = nameformat.format(snr.loc.galactic.l.value,snr.loc.galactic.b.value)
           snr.known = True
           snrs.append(snr)

    # Get the candidate SNRs from my text file
    if file_len("candidates.txt"):
        tempcat = np.loadtxt("candidates.txt")
# Makes it work for single-line files
        if len(tempcat.shape)==1:
           tempcat = [tempcat]

        for row in tempcat:
           snr = SNR()
           snr.loc = SkyCoord(ra = row[0]*u.degree, dec = row[1]*u.degree, frame='icrs')
           snr.maj = row[2]
           snr.min = row[2]
           snr.name  = nameformat.format(snr.loc.galactic.l.value,snr.loc.galactic.b.value)
           snr.known = False
           snrs.append(snr)

    print "Read "+str(len(snrs))+" SNRs from text files"
    return snrs


