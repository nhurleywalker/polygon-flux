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
        self.bkg = None # A dictionary of frequencies to background flux densities
        self.rms = None # A dictionary of frequencies to local RMS measurements
        self.nbeams = None # A dictionary of frequencies to number of PSFs fit over
        self.fit = None # A dictionary of "flux" -> total flux at 150 MHz; "alpha" -> fitted alpha; "chi2red" -> reduced chi2

def gaussian2d(x, y, mux, muy, sigmax, sigmay, theta):
   a = np.cos(theta)**2 / (2*sigmax**2) + np.sin(theta)**2 / (2*sigmay**2)
   b = -np.sin(2*theta) / (4*sigmax**2) + np.sin(2*theta) / (4*sigmay**2)
   c = np.sin(theta)**2 / (2*sigmax**2) + np.cos(theta)**2 / (2*sigmay**2)
   g = np.exp(-(a*(x-mux)**2 + 2*b*(x-mux)*(y-muy) + c*(y-muy)**2))
   return g

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0
# Read in the SNRs from files
# Get the known SNRs from Dave Green's "DAG.txt"
#                                          l b HH MM SS DD MM 
def read_snrs():
    if not os.path.exists("DAG.txt") and not os.path.exists("candidates.txt"):
        print "Neither DAG.txt nor candidates.txt exist so I don't know what to search for."
        sys.exit(1)
        
    snrs = []

    if is_non_zero_file("DAG.txt"):
        tempcat = np.genfromtxt("DAG.txt", delimiter=(6,6,4,3,3,5,3), comments="#")
        # Elliptical radii
        sizestrs = np.genfromtxt("DAG.txt",usecols=7,dtype=str, comments="#")
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
    if is_non_zero_file("candidates.txt"):
        tempcat = np.loadtxt("candidates.txt", comments="#")
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
    print fitsfile
    hdu = fits.open(fitsfile)

# Set any NaN areas to zero or the interpolation will fail
    hdu[0].data[np.isnan(hdu[0].data)] = 0.0

# Get vital stats of the fits file
    bmaj = hdu[0].header["BMAJ"]
    bmin = hdu[0].header["BMIN"]
    bpa = hdu[0].header["BPA"]
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
# Array of booleans of length len(indexes)
    inside = path.contains_points(indexes)

# New feature: nuke sources first, then do the rest
    for x,y in zip(local_sources.x,local_sources.y):
        if path.contains_points([[x,y]]):
            grid_x, grid_y = np.mgrid[0:xmax:1, 0:ymax:1]
#            radius = bmaj/pix2deg
# Exclude all indices inside the source from the interpolation
# TODO: maybe change this so it's more the shape of the PSF, to avoid wrecking perfectly useable bits of the image?
            print x, y, bmaj/(2*pix2deg), bmin/(2*pix2deg), np.radians(bpa)
            g2d = gaussian2d(grid_x, grid_y, y, x, bmaj/(2*pix2deg), bmin/(2*pix2deg), np.radians(180.-bpa))
            g = np.ravel(g2d)
            useforinterp = np.where(g < 0.1)
#            notforinterp = np.where(g > 0.1)
#            notforinterp = np.where(indexes.T[0]>(x-radius))
#            notforinterp = np.intersect1d(notforinterp, np.where(indexes.T[0]<(x+radius)))
#            notforinterp = np.intersect1d(notforinterp, np.where(indexes.T[1]<(y+radius)))
#            notforinterp = np.intersect1d(notforinterp, np.where(indexes.T[1]>(y-radius)))
#            useforinterp = np.setdiff1d(np.where(indexes),notforinterp)
            flux_values = []
            for ix in indexes[useforinterp]:
                flux_values.append(hdu[0].data[ix[1],ix[0]])
            interp = griddata(indexes[useforinterp], flux_values, (grid_x, grid_y), method="linear")
            hdu[0].data = interp.T
#            data_patch = hdu[0].data[np.floor(y-bmaj/pix2deg).astype(int):np.ceil(y+bmaj/pix2deg).astype(int),np.floor(x-bmaj/pix2deg).astype(int):np.ceil(x+bmaj/pix2deg).astype(int)]
#            peak_flux = np.max(data_patch) # In Jy/pix, but we assume the sources are unresolved
   # Save gaussian as a FITS file
#    header_new = hdu[0].header
#    new = fits.PrimaryHDU(g2d,header=header_new) #create new hdu
#    outname = fitsfile.replace(".fits","_gauss.fits")
#    newlist = fits.HDUList([new]) #create new hdulist
#    newlist.writeto(outname,overwrite=True)

   # Save source-interpolated as a FITS file
    header_new = hdu[0].header
        # Add a keyword so we know this one has been backgrounded
    header_new["INTERP"] = "TRUE"
    new = fits.PrimaryHDU(hdu[0].data,header=header_new) #create new hdu
    outname = fitsfile.replace(".fits","_interp.fits")
    newlist = fits.HDUList([new]) #create new hdulist
    newlist.writeto(outname,overwrite=True)

    if len(local_exclude.x):
        epath = mpltPath.Path(zip(local_exclude.x,local_exclude.y))
        ezone = epath.contains_points(indexes)

    # Total flux
    # Sadly astropy.io.fits is doing some sort of clever interpretation of this; have to do it in a loop instead
    # np.sum(hdu[0].data[indexes[np.where(inside)]])
    total_flux = 0.0
    for ix in indexes[np.where(inside)]:
       if not np.isnan(hdu[0].data[ix[1],ix[0]]):
           total_flux += hdu[0].data[ix[1],ix[0]] # In Jy/pix
#       else:
#           print "WARNING, NAN detected at x = {0}, y={1}!".format(ix[0], ix[1])
       
    total_flux = total_flux * (pix2deg**2) / beamvolume
    npix = indexes[np.where(inside)].shape[0]
    nbeams = npix * (pix2deg**2) / beamvolume

# Use a "ring" around the SNR for background, instead of interpolating and averaging over the interpolation
    # See https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.ndimage.morphology.binary_dilation.html
    # Note that we have to reshape the array back into 2D, otherwise it is just a long 1D list

    inner = np.reshape(ndimage.binary_dilation(np.reshape(inside,(xmax,ymax)),iterations=4),xmax*ymax)
    outer = np.reshape(ndimage.binary_dilation(np.reshape(inside,(xmax,ymax)),iterations=20),xmax*ymax)
# XOR operation: want to save the area where outer is True but inner is False,
# and exclude the areas where both are True (the SNR) and both are False (the rest of the image)
    bkger = [ a != b for a,b in zip(inner,outer)]
    bkg_indices = np.where(bkger)

# Exclude the areas the user flagged as unsuitable for a background fit
    if len(local_exclude.x):
        ex_indices = np.where(np.logical_not(ezone))
        useforbkg = np.intersect1d(bkg_indices, ex_indices)
    else:
        useforbkg = bkg_indices

   # Save mask as a FITS file
    header_new = hdu[0].header
        # Add a keyword so we know this one has been backgrounded
    header_new["MASK"] = "TRUE"
    mask = np.copy(hdu[0].data)
    mask[:,:] = np.nan
    # New feature: Now subtract the mean of the ring from the total_flux
    bkg_list = []
    for ix in indexes[useforbkg]:
       bkg_list.append(hdu[0].data[ix[1],ix[0]]) # In Jy/pix
       mask[ix[1],ix[0]] = hdu[0].data[ix[1],ix[0]] # Unset the NaNs of the mask
# Average background level in Jy/pix
    bkg_flux = np.mean(bkg_list)
# RMS of this area in Jy/pix
    rms = np.std(bkg_list)

    new = fits.PrimaryHDU(mask,header=header_new) #create new hdu
    outname = fitsfile.replace(".fits","_mask.fits")
    newlist = fits.HDUList([new]) #create new hdulist
    newlist.writeto(outname,overwrite=True)

#    source_flux = 0.0
#    for x,y in zip(local_sources.x,local_sources.y):
#       if path.contains_points([[x,y]]):
#    # Search for the local maximum within a beam-width
#    # Potential future improvement: modify this so it only does it for the white image, and then uses that location for all sources
#           data_patch = hdu[0].data[np.floor(y-bmaj/pix2deg).astype(int):np.ceil(y+bmaj/pix2deg).astype(int),np.floor(x-bmaj/pix2deg).astype(int):np.ceil(x+bmaj/pix2deg).astype(int)]
#           peak_flux = np.max(data_patch) # In Jy/pix, but we assume the sources are unresolved
##
##           bkg_index = np.unravel_index(data_patch.argmax(), data_patch.shape)
##           interp_patch = interp.T[np.floor(y-bmaj/pix2deg).astype(int):np.ceil(y+bmaj/pix2deg).astype(int),np.floor(x-bmaj/pix2deg).astype(int):np.ceil(x+bmaj/pix2deg).astype(int)]
##           bkg_at_src = interp_patch[bkg_index]
##           source_flux += (peak_flux - bkg_at_src)
#
## Subtract the average surface brightness of the remnant from the source flux
#           source_flux += (peak_flux - (total_flux / nbeams) )

    final_flux = total_flux - (bkg_flux * nbeams) # - source_flux

    #print "Number of beams searched: {0} \n Total flux density (Jy): {1} \n Total source flux density (Jy): {2} \n Background flux density (Jy): {3}\n Number of beams interpolated over after finding sources: {4}\n Final flux density (Jy): {5}".format(nbeams, total_flux, source_flux, bkg_flux, len(sources.x), final_flux)
    print "Number of beams searched: {0} \n Total flux density (Jy): {1} \n Background flux density (Jy): {2}\n Number of beams interpolated over after finding sources: {3}\n Final flux density (Jy): {4}".format(nbeams, total_flux, bkg_flux, len(sources.x), final_flux)

    #return final_flux, total_flux, source_flux, bkg_flux, rms
    return final_flux, total_flux, bkg_flux, rms, nbeams
