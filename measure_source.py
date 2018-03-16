#!/usr/bin/env python

__author__ = "Natasha Hurley-Walker"
__date__ = "16/03/2018"

import os
import sys
import shutil

import numpy as np
from astropy.io import fits
from astropy import wcs
from scipy.interpolate import griddata
from scipy import ndimage
# Need a least-squares estimator that gives a useable error estimate
from scipy.optimize import leastsq

import matplotlib
#matplotlib.use('Agg') # So does not use display -- only good if just making plots
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.path as mpltPath
from matplotlib.widgets import Slider, RadioButtons

import matplotlib.mlab as mlab
from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['serif']})

import argparse

import defs
from defs import Coords
from defs import SNR
from defs import read_snrs
from defs import PolyPick
from defs import create_index_array
from defs import find_fluxes

# Line markers and colors
firstsnr = { "marker" : "x" , "linestyle" : "None", "color" : "white" }
restsnr = { "marker" : "None" , "linestyle" : "-", "color" : "white" }
srcmark = { "marker" : "o" , "linestyle" : "None", "color" : "yellow" }
firstexc = { "marker" : "x" , "linestyle" : "None", "color" : "blue" }
restexc = { "marker" : "None" , "linestyle" : "--", "color" : "blue" }
reticsnr = { "marker" : "None" , "linestyle" : ":", "color" : "green" }

def poly_plot(fitsfile,makeplots):
    print "Fitting polygons to "+fitsfile
    hdu = fits.open(fitsfile)
    data = hdu[0].data
    xmax = hdu[0].header["NAXIS1"]
    ymax = hdu[0].header["NAXIS2"]
    w = wcs.WCS(hdu[0].header)
    try:
        pix2deg = hdu[0].header["CD2_2"]
    except KeyError:
        pix2deg = hdu[0].header["CDELT2"]

    def update(val):
        vmin = svmin.val
        vmax = svmax.val
        img.set_clim([svmin.val, svmax.val])
        fig.canvas.draw()
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([0.05,0.15,0.8,0.8])
    img = ax.imshow(data, cmap="gray")
    cbaxes = fig.add_axes([0.85, 0.1, 0.03, 0.85]) 
    cb = plt.colorbar(img, cax = cbaxes, orientation="vertical")  

    axcolor = 'white'
    axmin = fig.add_axes([0.05, 0.05, 0.2, 0.02], axisbg=axcolor)
    axmax  = fig.add_axes([0.44, 0.05, 0.2, 0.02], axisbg=axcolor)

    vmin0 = np.min(data)
    vmax0 = np.max(data)
    # include the slider: cribbed from https://stackoverflow.com/questions/5611805/using-matplotlib-slider-widget-to-change-clim-in-image
    svmin = Slider(axmin, "vmin", -5.0, 1.0, valinit=vmin0)
    svmax = Slider(axmax, "vmax", -1.0, 10.0, valinit=vmax0)

    svmin.on_changed(update)
    svmax.on_changed(update)

    fig.suptitle('Left-click = select source, middle-click = source subtract, right-click = select region to exclude from background calculation')

    # If you don't do this, it automatically zooms out when you first click on the plot
    ax.set_xlim(0,xmax)
    ax.set_ylim(0,ymax)
    polypick = PolyPick(ax)

    plt.show()

    # Export pixel co-ordinates as WCS co-ordinates in order to use on any arbitrary image, and save for later
    polygon = Coords()
    sources = Coords()
    exclude = Coords()
    if len(polypick.points.x):
        sources.x, sources.y = w.wcs_pix2world(zip(polypick.points.x,polypick.points.y),0).transpose()
    if len(polypick.exclude.x):
        exclude.x, exclude.y = w.wcs_pix2world(zip(polypick.exclude.x,polypick.exclude.y),0).transpose()

# Now I have the co-ordinates...
# Go measure the flux densities
    if len(polypick.coords.x):
        polygon.x, polygon.y = w.wcs_pix2world(zip(polypick.coords.x,polypick.coords.y),0).transpose()
        find_fluxes(polygon, sources, exclude, fitsfile)
# And make nice plots
        if makeplots is True:
            make_single_plot(polygon, sources, exclude, fitsfile)
    else:
        print "Failed to draw any points; exiting."

def make_single_plot(polygon, sources, exclude, fitsfile):
        print "Making attractive FITS image plot for "+fitsfile
# Load image data
        hdu = fits.open(fitsfile)
        hdu_interp = fits.open(fitsfile.replace(".fits","_interp.fits"))
        data = hdu[0].data
        data_interp = hdu_interp[0].data
## Get relevant header info
## Later: add bmaj bmin as ellipse
        xmax = hdu[0].header["NAXIS1"]
        ymax = hdu[0].header["NAXIS2"]
## Transform polygons into local pixel co-ordinates
        w = wcs.WCS(hdu[0].header)
        local_polygon = Coords()
        local_polygon.x, local_polygon.y = w.wcs_world2pix(zip(polygon.x,polygon.y),0).transpose()
        local_sources = Coords()
        if len(sources.x):
            local_sources.x, local_sources.y = w.wcs_world2pix(zip(sources.x,sources.y),0).transpose()
        local_exclude = Coords()
        if len(exclude.x):
            local_exclude.x, local_exclude.y = w.wcs_world2pix(zip(exclude.x,exclude.y),0).transpose()

# Using http://docs.astropy.org/en/stable/visualization/wcsaxes/index.html
        fig = plt.figure()
        ax = fig.add_subplot(121, projection=w)
        ax_interp = fig.add_subplot(122, projection=w)
        img = ax.imshow(data, cmap="gray",origin="lower")
        img_interp = ax_interp.imshow(data_interp, cmap="gray",origin="lower")
        for a in ax, ax_interp:
            overlay = a.get_coords_overlay("galactic")
            overlay.grid(axes = a, color='white', ls='dotted')
            a.set_xlim(0,xmax)
            a.set_ylim(0,ymax)
            if len(local_polygon.x):
                a.plot(local_polygon.x,local_polygon.y,**restsnr)
            if len(local_sources.x):
                a.plot(local_sources.x,local_sources.y,**srcmark)
            if len(local_exclude.x):
                a.plot(local_exclude.x,local_exclude.y,**restexc)
        fig.savefig(fitsfile.replace(".fits",".png"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("files to read")
    group1.add_argument('--fitsfile', dest='fitsfile', default=None, \
                        help="Low-resolution file to read")
    group2 = parser.add_argument_group("functions to perform")
    group2.add_argument('--plot', dest='makeplots', default=False, action='store_true', \
                        help="Make png plots of the FITS image (default = False)")

    options = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit()

    # Perform the fitting
    poly_plot(options.fitsfile,options.makeplots)
