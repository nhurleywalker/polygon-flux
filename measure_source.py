#!/usr/bin/env python

__author__ = "Natasha Hurley-Walker"
__date__ = "16/03/2018"

import os
import sys
import shutil

import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.visualization import PercentileInterval

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

def poly_plot(fitsfile, makeplots, bkgtype="mean", output=False):
    print("Fitting polygons to "+fitsfile)
    hdu = fits.open(fitsfile)
# Remove degenerate axes
    data = np.squeeze(np.squeeze(hdu[0].data))
    xmax = hdu[0].header["NAXIS1"]
    ymax = hdu[0].header["NAXIS2"]
    w = wcs.WCS(hdu[0].header, naxis=2)
    try:
        pix2deg = hdu[0].header["CD2_2"]
    except KeyError:
        pix2deg = hdu[0].header["CDELT2"]

# Percentile interval for image normalisation
    pct = 99.0
    interval = PercentileInterval(pct)
    vmin, vmax = interval.get_limits(data)

    def update(val):
#        vmin = svmin.val
#        vmax = svmax.val
        img.set_clim([svmin.val, svmax.val])
        fig.canvas.draw()
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([0.05,0.15,0.8,0.8])
    img = ax.imshow(data, cmap="gray", vmin=vmin, vmax=vmax)
    cbaxes = fig.add_axes([0.85, 0.1, 0.03, 0.85]) 
    cb = plt.colorbar(img, cax = cbaxes, orientation="vertical")  

    axcolor = 'white'
    try:
        axmin = fig.add_axes([0.05, 0.05, 0.2, 0.02], axisbg=axcolor)
        axmax  = fig.add_axes([0.44, 0.05, 0.2, 0.02], axisbg=axcolor)
    except AttributeError:
        axmin = fig.add_axes([0.05, 0.05, 0.2, 0.02], facecolor=axcolor)
        axmax  = fig.add_axes([0.44, 0.05, 0.2, 0.02], facecolor=axcolor)

    vmin0 = np.min(data)
    vmax0 = np.max(data)
    # include the slider: cribbed from https://stackoverflow.com/questions/5611805/using-matplotlib-slider-widget-to-change-clim-in-image
    svmin = Slider(axmin, "vmin", vmin, vmax, valinit=vmin)
    svmax = Slider(axmax, "vmax", vmin, vmax, valinit=vmax)

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
#        sources.x, sources.y = w.wcs_pix2world(zip(polypick.points.x,polypick.points.y),0).transpose()
        sources.x, sources.y = w.wcs_pix2world(polypick.points.x,polypick.points.y,0)
    if len(polypick.exclude.x):
#        exclude.x, exclude.y = w.wcs_pix2world(zip(polypick.exclude.x,polypick.exclude.y),0).transpose()
        exclude.x, exclude.y = w.wcs_pix2world(polypick.exclude.x,polypick.exclude.y,0)

# Now I have the co-ordinates...
# Go measure the flux densities
    if len(polypick.coords.x):
        polygon.x, polygon.y = w.wcs_pix2world(polypick.coords.x,polypick.coords.y,0)
        find_fluxes(polygon, sources, exclude, fitsfile, bkgtype, output)

# And make nice plots
        if makeplots is True:
            make_single_plot(polygon, sources, exclude, fitsfile)
    else:
        print("Failed to draw any points; exiting.")

def make_single_plot(polygon, sources, exclude, fitsfile):
        print("Making attractive FITS image plot for "+fitsfile)
# Load image data
        hdu = fits.open(fitsfile)
        hdu_mask = fits.open(fitsfile.replace(".fits","_mask.fits"))
        data = np.squeeze(np.squeeze(hdu[0].data))
        data_mask = hdu_mask[0].data
## Get relevant header info
## Later: add bmaj bmin as ellipse
        xmax = hdu[0].header["NAXIS1"]
        ymax = hdu[0].header["NAXIS2"]
## Transform polygons into local pixel co-ordinates
        w = wcs.WCS(hdu[0].header, naxis=2)
        local_polygon = Coords()
        local_polygon.x, local_polygon.y = w.wcs_world2pix(polygon.x,polygon.y,0)
        local_sources = Coords()
        if len(sources.x):
            local_sources.x, local_sources.y = w.wcs_world2pix(sources.x,sources.y,0)
        local_exclude = Coords()
        if len(exclude.x):
            local_exclude.x, local_exclude.y = w.wcs_world2pix(exclude.x,exclude.y,0)

# Make image visible by using interval normalisation
# Percentile interval for image normalisation
        pct = 98.5
        interval = PercentileInterval(pct)
        vmin, vmax = interval.get_limits(data)

# Using http://docs.astropy.org/en/stable/visualization/wcsaxes/index.html
        fig = plt.figure()
# Unadorned panel
        ax = fig.add_subplot(121, projection=w)
# Panel with the backgrounding etc shown
        ax_mask = fig.add_subplot(122, projection=w)
        img_l = ax.imshow(data, cmap="gray",origin="lower", vmin=vmin, vmax=vmax)
        img_r = ax_mask.imshow(data, cmap="gray",origin="lower", vmin=vmin, vmax=vmax)
        img_mask = ax_mask.imshow(data_mask, cmap="cool_r",origin="lower", vmin=9999., vmax=9999., alpha=0.6)
#        img = ax.imshow(data, cmap="gray",origin="lower", norm = matplotlib.colors.LogNorm(vmin=0.0, vmax=np.nanmax(data)))
#        img_mask = ax_mask.imshow(data_mask, cmap="gray",origin="lower", norm = matplotlib.colors.LogNorm(vmin=0.0, vmax=np.nanmax(data)))
#        img = ax.imshow(data, cmap="gray",origin="lower", norm = matplotlib.colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=np.nanmin(data), vmax=np.nanmax(data)))
#        img_mask = ax_mask.imshow(data_mask, cmap="gray",origin="lower", norm = matplotlib.colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=np.nanmin(data), vmax=np.nanmax(data)))
        for a in ax, ax_mask:
            overlay = a.get_coords_overlay("galactic")
            overlay.grid(axes = a, color='white', ls='dotted')
            a.set_xlim(0,xmax)
            a.set_ylim(0,ymax)
            lon = a.coords['ra']
            lon.set_axislabel("Right Ascension (J2000)")
            lon.set_major_formatter('hh:mm')
# Nice Dec labels for left panel
        lat = ax.coords['dec']
        lat.set_axislabel("Declination (J2000)")
        lat.set_major_formatter('dd:mm')
# Overlay plot of excluded regions and selected region
        if len(local_polygon.x):
            ax_mask.plot(local_polygon.x,local_polygon.y,**restsnr)
        if len(local_sources.x):
            ax_mask.plot(local_sources.x,local_sources.y,**srcmark)
        if len(local_exclude.x):
            ax_mask.plot(local_exclude.x,local_exclude.y,**restexc)
# No Dec labels for right panel
        lat = ax_mask.coords['dec']
        lat.set_ticklabel_visible(False)
        lat.set_ticks_visible(False)
        fig.savefig(fitsfile.replace(".fits",".png"), bbox_inches="tight")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("files to read")
    group1.add_argument('--fitsfile', dest='fitsfile', default=None, \
                        help="Low-resolution file to read")
    group2 = parser.add_argument_group("Generate plots and images")
    group2.add_argument('--plot', dest='makeplots', default=False, action='store_true', \
                        help="Make png plots of the FITS image (default = False)")
    group2.add_argument('--output', dest='output', default=False, action='store_true', \
                        help="Make FITS images of the background-subtracted image etc (default = False)")
    group3 = parser.add_argument_group("Options")
    group2.add_argument('--bkgtype', dest='bkgtype', default="mean", \
                        help="Use mean background ('mean'), or an interpolated 2D plane ('interp') (default = 'mean')")


    options = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit()

# Test that the image is square -- otherwise the dilation to find the border fails
    square = True
    header = fits.getheader(options.fitsfile)
    if header["NAXIS"] == 4:
        if header["NAXIS3"] != header["NAXIS4"]:
            square = False
    else:
        if header["NAXIS1"] != header["NAXIS2"]:
            square = False

    if square is True:
        # Perform the fitting
        poly_plot(options.fitsfile,options.makeplots, bkgtype=options.bkgtype, output=options.output)
    else:
        print("Please provide a square FITS image.")
