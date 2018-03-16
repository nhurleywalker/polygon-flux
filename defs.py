#!/usr/bin/env python

__author__ = "Natasha Hurley-Walker"
__date__ = "10/03/2018"

import os
import sys
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u

from scipy.interpolate import griddata
from scipy import ndimage

import matplotlib
from matplotlib import pyplot as plt
import matplotlib.path as mpltPath

import numpy as np

from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['serif']})

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

def create_index_array(hdu):
    # Initialise an array of co-ordinate pairs for the full array
    indexes = np.empty(((hdu[0].data.shape[0])*(hdu[0].data.shape[1]),2),dtype=int)
    idx = np.array([ (j,0) for j in xrange(hdu[0].data.shape[1])])
    j=hdu[0].data.shape[1]
    for i in xrange(hdu[0].data.shape[0]):
        idx[:,1]=i
        indexes[i*j:(i+1)*j] = idx
    return indexes

# Hacked from https://matplotlib.org/users/event_handling.html
class PolyPick:
    def __init__(self, ax=None):
        if ax is None:
            self.ax = plt.gca()
        else:
            self.ax = ax
        self.cid = ax.figure.canvas.mpl_connect('button_press_event', self)
# Lines (region to include in remnant)
        self.coords = Coords()
# Points (to find and exclude point sources)
        self.points = Coords()
# Exclude (region to exclude from background fit
        self.exclude = Coords()

    def __call__(self, event):
        if event.inaxes!=self.ax.axes: return
        if event.button == 1:
            print "Selected coords {0:5.4f}, {1:5.4f}".format(event.xdata,event.ydata)
            self.coords.x.append(event.xdata)
            self.coords.y.append(event.ydata)
        elif event.button == 2:
            print "Selecting source at {0:5.4f}, {1:5.4f} for fitting".format(event.xdata,event.ydata)
            self.points.x.append(event.xdata)
            self.points.y.append(event.ydata)
        elif event.button == 3:
            print "Selecting coords at {0:5.4f}, {1:5.4f} to remove from background fit".format(event.xdata,event.ydata)
            self.exclude.x.append(event.xdata)
            self.exclude.y.append(event.ydata)

# Remnant region
# First line: draw a "x" to show it's the first one
        if len(self.coords.x) == 1:
            line, = self.ax.plot(self.coords.x, self.coords.y,**firstsnr)
# Second or further line: draw the whole line
        else:
            line, = self.ax.plot(self.coords.x, self.coords.y,**restsnr)
# Sources
        line, = self.ax.plot(self.points.x, self.points.y,**srcmark)
# Exclusion zone
        if len(self.exclude.x) == 1:
            line, = self.ax.plot(self.exclude.x, self.exclude.y,**firstexc)
# Second or further line: draw the whole line
        else:
            line, = self.ax.plot(self.exclude.x, self.exclude.y,**restexc)

# This makes it plot without needing to change focus back to the terminal
        self.ax.figure.canvas.draw()
#    def sources():
#         return self.points.x, self.points.y
#    def coords():
#         return self.coords.x, self.coords.y

def find_fluxes(polygon, sources, exclude, fitsfile):#, export):
    hdu = fits.open(fitsfile)
    # Transform the polygon and source arrays into local pixel co-ordinates for further operations
    w = wcs.WCS(hdu[0].header)
    local_polygon = Coords()
    local_polygon.x, local_polygon.y = w.wcs_world2pix(zip(polygon.x,polygon.y),0).transpose()
    local_sources = Coords()
    if len(sources.x):
        local_sources.x, local_sources.y = w.wcs_world2pix(zip(sources.x,sources.y),0).transpose()
    local_exclude = Coords()
    if len(exclude.x):
        local_exclude.x, local_exclude.y = w.wcs_world2pix(zip(exclude.x,exclude.y),0).transpose()
    
    indexes = create_index_array(hdu)
    # Adapted from https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python
    path = mpltPath.Path(zip(local_polygon.x,local_polygon.y))
    inside = path.contains_points(indexes)

    if len(local_exclude.x):
        epath = mpltPath.Path(zip(local_exclude.x,local_exclude.y))
        ezone = epath.contains_points(indexes)

    # Total flux
    # Sadly astropy.io.fits is doing some sort of clever interpretation of this; have to do it in a loop instead
    # np.sum(hdu[0].data[indexes[np.where(inside)]])
    total_flux = 0.0
    for ix in indexes[np.where(inside)]:
       total_flux += hdu[0].data[ix[1],ix[0]] # In Jy/pix
       
    bmaj = hdu[0].header["BMAJ"]
    bmin = hdu[0].header["BMIN"]
    xmax = hdu[0].header["NAXIS1"]
    ymax = hdu[0].header["NAXIS2"]
    try:
        pix2deg = hdu[0].header["CD2_2"]
    except KeyError:
        pix2deg = hdu[0].header["CDELT2"]
# Montaged images use PC instead of CD
    if pix2deg == 1.0:
        pix2deg = hdu[0].header["PC2_2"]
    beamvolume = (1.1331 * bmaj * bmin) # gaussian beam conversion

    total_flux = total_flux * (pix2deg**2) / beamvolume
    npix = indexes[np.where(inside)].shape[0]
    nbeams = npix * (pix2deg**2) / beamvolume

    # Feed the interpolation routine the slightly larger area, so we avoid interpolating directly on my SNR
    # See https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.ndimage.morphology.binary_dilation.html
    # Note that we have to reshape the array back into 2D, otherwise it is just a long 1D list
    # TODO: now that I'm using interpolated images, it might be more suitable to make the number of iterations dependent on the frequency in the fits header, because the red images are highly oversampled. But let's see if it makes any difference first.
    grown = np.reshape(ndimage.binary_dilation(np.reshape(inside,(xmax,ymax)),iterations=4),xmax*ymax)

    # Exclude the areas unsuitable for a background fit
    bkg_indices = np.where(np.logical_not(grown))
    if len(local_exclude.x):
        ex_indices = np.where(np.logical_not(ezone))
        useforinterp = np.intersect1d(bkg_indices, ex_indices)
    else:
        useforinterp = bkg_indices

    grid_x, grid_y = np.mgrid[0:xmax-1:1, 0:ymax-1:1]

    flux_values = []
#    for ix in indexes[np.where(np.logical_not(grown))]:
    for ix in indexes[useforinterp]:
       flux_values.append(hdu[0].data[ix[1],ix[0]])

    interp = griddata(indexes[useforinterp], flux_values, (grid_x, grid_y), method="linear")

    # TODO: restore the data values for the exclusion zone so that the interpolated fits file still has that area ready for plotting later

   # Save as a FITS file
    header_new = hdu[0].header
        # Add a keyword so we know this one has been backgrounded
    header_new["INTERP"] = "TRUE"
    new = fits.PrimaryHDU(interp.T,header=header_new) #create new hdu
    outname = fitsfile.replace(".fits","_interp.fits")
    newlist = fits.HDUList([new]) #create new hdulist
    newlist.writeto(outname,overwrite=True)

    # Now subtract what we interpolated from the total_flux
    bkg_flux = 0.0
    for ix in indexes[np.where(inside)]:
       bkg_flux += interp[ix[0],ix[1]] # In Jy/pix

    bkg_flux = bkg_flux * (pix2deg**2) / beamvolume

    source_flux = 0.0
    for x,y in zip(local_sources.x,local_sources.y):
       if path.contains_points([[x,y]]):
    # Search for the local maximum within a beam-width
           data_patch = hdu[0].data[np.floor(y-bmaj/pix2deg).astype(int):np.ceil(y+bmaj/pix2deg).astype(int),np.floor(x-bmaj/pix2deg).astype(int):np.ceil(x+bmaj/pix2deg).astype(int)]
           peak_flux = np.max(data_patch) # In Jy/pix, but we assume the sources are unresolved
    # Need the background level at this point
           bkg_index = np.unravel_index(data_patch.argmax(), data_patch.shape)
           interp_patch = interp.T[np.floor(y-bmaj/pix2deg).astype(int):np.ceil(y+bmaj/pix2deg).astype(int),np.floor(x-bmaj/pix2deg).astype(int):np.ceil(x+bmaj/pix2deg).astype(int)]
           bkg_at_src = interp_patch[bkg_index]
           source_flux += (peak_flux - bkg_at_src)
   # TODO: modify the SNR source location to match the peak flux location -- will need to modify the SNR object at the end of this function

    final_flux = total_flux - source_flux - bkg_flux

    print "Number of beams searched: {0} \n Total flux density (Jy): {1} \n Total source flux density (Jy): {2} \n Background flux density (Jy): {3}\n Number of beams masked after finding sources: {4}\n Final flux density (Jy): {5}".format(nbeams, total_flux, source_flux, bkg_flux, len(sources.x), final_flux)

    return final_flux
