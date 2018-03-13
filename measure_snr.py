#!/usr/bin/env python

__author__ = "Natasha Hurley-Walker"
__date__ = "10/03/2018"

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

# Line markers and colors
firstsnr = { "marker" : "x" , "linestyle" : "None", "color" : "white" }
restsnr = { "marker" : "None" , "linestyle" : "-", "color" : "white" }
srcmark = { "marker" : "o" , "linestyle" : "None", "color" : "yellow" }
firstexc = { "marker" : "x" , "linestyle" : "None", "color" : "blue" }
restexc = { "marker" : "None" , "linestyle" : "--", "color" : "blue" }
reticsnr = { "marker" : "None" , "linestyle" : ":", "color" : "green" }

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

def create_index_array(hdu):
    # Initialise an array of co-ordinate pairs for the full array
    indexes = np.empty(((hdu[0].data.shape[0])*(hdu[0].data.shape[1]),2),dtype=int)
    idx = np.array([ (j,0) for j in xrange(hdu[0].data.shape[1])])
    j=hdu[0].data.shape[1]
    for i in xrange(hdu[0].data.shape[0]):
        idx[:,1]=i
        indexes[i*j:(i+1)*j] = idx
    return indexes

# define our (line) fitting function
fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

# http://scipy-cookbook.readthedocs.org/items/FittingData.html
# Define function for calculating a power law
def powerlaw(x,amp,index):
    return amp * (x**index)

def fit_spectrum(freq_array,flux_array,flux_errors): #,plot):
    pinit = [3.0, -0.7]
    fit = leastsq(errfunc, pinit, args=(freq_array, flux_array, flux_errors), full_output=1)
    covar = fit[1]
    if covar is not None:
        P = fit[0]
        residual = errfunc(P,freq_array, flux_array, flux_errors)
        chi2red = sum(np.power(residual,2))/(len(freq_array)-len(pinit))
        amp = P[0]
        alpha = P[1]
    # Errors
        err_amp = np.sqrt(covar[0][0])
        err_alpha = np.sqrt(covar[1][1])
    else:
        chi2red=None
        alpha=None
        amp=None
        err_alpha=None
        err_amp=None
    return alpha, err_alpha, amp, err_amp, chi2red

# Allow me to display the figure in the centre of the screen, more pleasant than at a random location
def newfigure(num=None,**args):
       hfig = plt.figure(num,**args)
# TODO: Make this work for typical window managers
# Or fail gracefully in the case of not knowing which function to call
       plt.get_current_fig_manager().window.set_position("GTK_WIN_POS_CENTER")
       plt.get_current_fig_manager().window.resize(1000,500)
       return hfig

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

def poly_plots(snrs,export):
    snrs2 = []
    for snr in snrs:
        print "Fitting polygons and sources to "+snr.name
        white = fits.open("white/"+snr.name+".fits")
        red = fits.open("red/rpj/"+snr.name+".fits")
        green = fits.open("green/rpj/"+snr.name+".fits")
        blue = fits.open("blue/rpj/"+snr.name+".fits")
        white_data = white[0].data
        red_data = red[0].data
        green_data = green[0].data
        blue_data = blue[0].data
        rgb = np.dstack([red_data,green_data,blue_data])
        xmax = white[0].header["NAXIS1"]
        ymax = white[0].header["NAXIS2"]
        w = wcs.WCS(white[0].header)
        try:
            pix2deg = white[0].header["CD2_2"]
        except KeyError:
            pix2deg = white[0].header["CDELT2"]
        # Subtract minimum value (usually negative) in order to make the plotting nice
        rgb -= np.nanmin(rgb)
        # Divide by maximum value to normalise
        rgb /= np.nanmax(rgb)

        def update(val):
            vmin = svmin.val
            vmax = svmax.val
            img_white.set_clim([svmin.val, svmax.val])
            fig.canvas.draw()
# TODO: Turn this on or off
#        fig = plt.figure()
        fig = newfigure(figsize=(10,5))
        ax_white = fig.add_axes([0.05,0.15,0.4,0.8])
        ax_rgb = fig.add_axes([0.55,0.15,0.4,0.8])
        img_rgb = ax_rgb.imshow(rgb)
        img_white = ax_white.imshow(white_data, cmap="gray")
        cbaxes = fig.add_axes([0.45, 0.1, 0.03, 0.85]) 
        cb = plt.colorbar(img_white, cax = cbaxes, orientation="vertical")  

        axcolor = 'white'
        axmin = fig.add_axes([0.05, 0.05, 0.1, 0.02], axisbg=axcolor)
        axmax  = fig.add_axes([0.24, 0.05, 0.1, 0.02], axisbg=axcolor)

        vmin0 = np.min(white_data)
        vmax0 = np.max(white_data)
        # include the slider: cribbed from https://stackoverflow.com/questions/5611805/using-matplotlib-slider-widget-to-change-clim-in-image
        svmin = Slider(axmin, "vmin", -5.0, 1.0, valinit=vmin0)
        svmax = Slider(axmax, "vmax", -1.0, 10.0, valinit=vmax0)

        svmin.on_changed(update)
        svmax.on_changed(update)

        fig.suptitle('Left-click = select source, middle-click = source subtract, right-click = select region to exclude from background calculation')

        centerx, centery = w.wcs_world2pix(snr.loc.galactic.l.value,snr.loc.galactic.b.value,0)

        major, minor = snr.maj/(2*pix2deg) , snr.min/(2*pix2deg)
        for ax in ax_white, ax_rgb:
        # If you don't do this, it automatically zooms out when you first click on the plot
            ax.set_xlim(0,xmax)
            ax.set_ylim(0,ymax)
        # Plot rough size of SNR -- DAG doesn't include pa information so neither can I
            ax.plot([centerx,centerx],[centery-major,centery+major],**reticsnr)
            ax.plot([centerx-minor,centerx+minor],[centery,centery],**reticsnr)
        polypick = PolyPick(ax_white)

        rax = fig.add_axes([0.6, 0.0, 0.15, 0.10])
        radio = RadioButtons(rax, ("white", "rgb"), active=0)

        def changeax(axl):
            if axl == "white":
                polypick.ax = ax_white
                plt.sca(ax_white)
            elif axl == "rgb":
                polypick.ax = ax_rgb
                plt.sca(ax_rgb)
            fig.canvas.draw_idle()
        radio.on_clicked(changeax)

        plt.show()

        # Export pixel co-ordinates as WCS co-ordinates in order to use on any arbitrary image, and save for later
        polygon = Coords()
        sources = Coords()
        exclude = Coords()
        if len(polypick.coords.x):
            polygon.x, polygon.y = w.wcs_pix2world(zip(polypick.coords.x,polypick.coords.y),0).transpose()
        else:
# User has got bored of picking polygons and wants to exit this loop
            print "No polygon detected; leaving polygon-drawing mode."
            snrs = snrs2
            break
            
        if len(polypick.points.x):
            sources.x, sources.y = w.wcs_pix2world(zip(polypick.points.x,polypick.points.y),0).transpose()
        if len(polypick.exclude.x):
            exclude.x, exclude.y = w.wcs_pix2world(zip(polypick.exclude.x,polypick.exclude.y),0).transpose()
# TODO: fix this so it is in RA and Dec, not galactic co-ordinates!!

        snr.polygon = polygon
        snr.sources = sources
        snr.exclude = exclude
        snrs2.append(snr)

    if export:
        export_snrs(snrs)

    return snrs2 # Is this necessary?

def update_text_file(snr, outtype):
    if outtype == "latex":
        outformat = "{0} & {1:02.0f}:{2:02.0f}:{3:02.0f} & {4:+02.0f}:{5:02.0f}:{6:02.0f} & ${7:3.1f}\\times{8:3.1f}$ & ${9:3.2f}\pm{10:3.2f}$ & ${11:3.2f}\pm{12:3.2f}$ \\\\ \n"
        outvars = [snr.name, \
                snr.loc.ra.hms.h, snr.loc.ra.hms.m, snr.loc.ra.hms.s, \
                snr.loc.dec.dms.d, abs(snr.loc.dec.dms.m), abs(snr.loc.dec.dms.s), \
                60*snr.maj, 60*snr.min, \
                snr.fit["flux"], snr.fit["fluxerr"], snr.fit["alpha"], snr.fit["alphaerr"]]
        outputfile = "SNR_flux_densities.tex"
        if not os.path.exists(outputfile):
            with open(outputfile, 'w') as output_file:
               output_file.write("Name & RA & Dec & Size & Flux & $\\alpha$ \\\\ \n")
    elif outtype == "csv":
        outformat = "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9}\n"
        outvars = [snr.name, \
                snr.loc.ra.value, snr.loc.dec.value, \
                60*snr.maj, 60*snr.min, \
                snr.fit["flux"], snr.fit["fluxerr"], snr.fit["alpha"], snr.fit["alphaerr"], snr.fit["chi2red"]]
        outputfile = "SNR_flux_densities.csv"
        if not os.path.exists(outputfile):
            with open(outputfile, 'w') as output_file:
               output_file.write("#Name,RA,Dec,Major,Minor,Flux,Err_flux,Alpha,Err_Alpha\n")
# If the SNR has an entry, overwrite it
    written = False
    with open(outputfile, 'r') as input_file, open('new_file.txt', 'w') as output_file:
        for line in input_file:
            if snr.name in line:
                output_file.write(outformat.format(*outvars))
                written = True
            else:
                output_file.write(line)
    shutil.copy("new_file.txt", outputfile)
# If the SNR does not have an entry, append it to the file
    if written is False:
        with open(outputfile, 'a') as output_file:
            output_file.write(outformat.format(*outvars))
# TODO: sort the output file by some index -- galactic longitude then latitude?

def fit_fluxes(snrs, export):#, makeplots):
    normfreq = 150000000. # 150 MHz
    for snr in snrs:
        print "Fitting spectrum to flux densities from "+snr.name
        freqs = []
        fluxes = []
        # Loop over the colors
        colors = ["white", "red", "green", "blue"]
        for color in colors:
            if color == "white":
                fitsfile = color+"/"+snr.name+".fits"
            else:
                fitsfile = color+"/rpj/"+snr.name+".fits"
            hdu = fits.open(fitsfile)
            freqs.append(hdu[0].header["FREQ"])
            fluxes.append(find_fluxes(snr.polygon, snr.sources, snr.exclude, fitsfile))#, export))
# Save to SNR object for later
        snr.flux = dict(zip(freqs, fluxes))

    # Obviously get this from the header in future
#        freqs = [88, 118, 154, 200]
        logfreqs = np.log([f/normfreq for f in freqs])
        logfluxes = np.log(fluxes)
        
#        print np.log(freqs)
#        print np.log(fluxes)
# 10% errors to start with
        fluxerrors = [0.1*s for s in fluxes]
        alpha, err_alpha, amp, err_amp, chi2red = fit_spectrum(logfreqs,logfluxes,0.1)

        # frequencies for plotting
        plotfreqs = np.linspace(0.9*np.min(freqs),1.1*np.max(freqs),200)

        if alpha:
            flux150 = np.exp(amp)
            flux150err = err_amp * flux150
            print "Spectrum fit, making plot for "+snr.name
            outpng = ("./spectra/"+snr.name+".png")
            # Plot
            example = plt.figure(figsize=(10,5))
            # LaTeX-safe
            example.suptitle(snr.name.replace("_","\,"))
            ax1 = example.add_subplot(1,2,1)
            ax1.plot(plotfreqs, powerlaw(plotfreqs, flux150/(normfreq**alpha), alpha),label="alpha=${0:3.1f}$ ".format(alpha))     # Scipy Fit
            ax1.errorbar(freqs, fluxes, yerr=fluxerrors, fmt='k.')  # Data
            ax1.set_ylabel("S / Jy")
            ax1.set_xlabel("Frequency / MHz")
            ax1.legend()
            ax2=example.add_subplot(1,2,2)
            ax2.loglog(plotfreqs, powerlaw(plotfreqs, flux150/(normfreq**alpha), alpha),label="alpha=${0:3.1f}$ ".format(alpha))     # Scipy Fit
            ax2.errorbar(freqs, fluxes, yerr=fluxerrors, fmt='k.')  # Data
    #        ax2.set_xlim(left=xmin01,right=xmax01)
    #       ax2.set_ylim([min(np.exp(flux_array)),max(np.exp(flux_array))])
            ax2.set_ylabel("S / Jy")
            ax2.set_xlabel("Frequency / MHz")
            ax2.legend()
            example.savefig(outpng)
            example.clf()
        
            snr.fit = {"flux" : flux150, "fluxerr" : flux150err, "alpha" : alpha, "alphaerr" : err_alpha, "chi2red" : chi2red}
# Add table output
            for outtype in ["latex", "csv"]:
                update_text_file(snr, outtype)
    return snrs

def export_snrs(snrs):
    print "Exporting SNR catalogue"
    output_file = "SNR_catalogue.npz"
    np.savez(output_file,snrs=snrs)
    # Also need to add a pretty output format table for reference, maybe a FITS or VO?

def import_snrs():
    print "Importing SNR catalogue"
    input_file = "SNR_catalogue.npz"
    if os.path.exists(input_file):
        data = np.load(input_file)
        snrs = data["snrs"].tolist()
    else:
        print "Cannot find existing SNR_catalogue.npz!"
        snrs = []
    print("Imported {0} SNRs".format(len(snrs)))
    return snrs

def make_plots(snrs):
    for snr in snrs:
        print "Making attractive FITS image plot for "+snr.name
# Load image data
        white = fits.open("white/"+snr.name+".fits")
        white_interp = fits.open("white/"+snr.name+"_interp.fits")
        red = fits.open("red/rpj/"+snr.name+".fits")
        red_interp = fits.open("red/rpj/"+snr.name+"_interp.fits")
        green = fits.open("green/rpj/"+snr.name+".fits")
        green_interp = fits.open("green/rpj/"+snr.name+"_interp.fits")
        blue = fits.open("blue/rpj/"+snr.name+".fits")
        blue_interp = fits.open("blue/rpj/"+snr.name+"_interp.fits")
        white_data = white[0].data
        white_interp_data = white_interp[0].data
        red_data = red[0].data
        red_interp_data = red_interp[0].data
        green_data = green[0].data
        green_interp_data = green_interp[0].data
        blue_data = blue[0].data
        blue_interp_data = blue_interp[0].data
        rgb = np.dstack([red_data,green_data,blue_data])
        rgb_interp = np.dstack([red_interp_data,green_interp_data,blue_interp_data])
        # Subtract minimum value (usually negative) in order to make the plotting nice
        rgb -= np.min([np.nanmin(rgb),np.nanmin(rgb_interp)])
        rgb_interp -= np.min([np.nanmin(rgb),np.nanmin(rgb_interp)])
        # Divide by maximum value to normalise
        rgb /= np.max([np.nanmax(rgb),np.nanmax(rgb_interp)])
        rgb_interp /= np.max([np.nanmax(rgb),np.nanmax(rgb_interp)])
## Get relevant header info
## Later: add bmaj bmin as ellipse
        xmax = white[0].header["NAXIS1"]
        ymax = white[0].header["NAXIS2"]
## Transform snr polygons into local pixel co-ordinates
## All colours should be on the same projection thanks to Montage so we only have to do this once
        w = wcs.WCS(white[0].header)
        local_polygon = Coords()
        if len(snr.polygon.x):
            local_polygon.x, local_polygon.y = w.wcs_world2pix(zip(snr.polygon.x,snr.polygon.y),0).transpose()
        local_sources = Coords()
        if len(snr.sources.x):
            local_sources.x, local_sources.y = w.wcs_world2pix(zip(snr.sources.x,snr.sources.y),0).transpose()
        local_exclude = Coords()
        if len(snr.exclude.x):
            local_exclude.x, local_exclude.y = w.wcs_world2pix(zip(snr.exclude.x,snr.exclude.y),0).transpose()
        print local_exclude.x, local_exclude.y
        print local_sources.x, local_sources.y

# Using http://docs.astropy.org/en/stable/visualization/wcsaxes/index.html
        fig = plt.figure()
        ax_white = fig.add_subplot(221, projection=w)
        ax_white_interp = fig.add_subplot(222, projection=w)
        ax_rgb = fig.add_subplot(223, projection=w)
        ax_rgb_interp = fig.add_subplot(224, projection=w)
        img_white = ax_white.imshow(white_data, cmap="gray",origin="lower")
        img_white_interp = ax_white_interp.imshow(white_interp_data, cmap="gray",origin="lower")
        img_rgb = ax_rgb.imshow(rgb,origin="lower")
        img_rgb_interp = ax_rgb_interp.imshow(rgb_interp,origin="lower")
        for ax in ax_white, ax_white_interp, ax_rgb, ax_rgb_interp:
            overlay = ax.get_coords_overlay("galactic")
            overlay.grid(axes = ax, color='white', ls='dotted')
            ax.set_xlim(0,xmax)
            ax.set_ylim(0,ymax)
            if len(local_polygon.x):
                ax.plot(local_polygon.x,local_polygon.y,**restsnr)
            if len(local_sources.x):
                ax.plot(local_sources.x,local_sources.y,**srcmark)
            if len(local_exclude.x):
                ax.plot(local_exclude.x,local_exclude.y,**restexc)
        fig.savefig("plots/"+snr.name+".png")

#       figure out what plots I actually want:
#       some combination of:
#             - white vs rgb
#             - with annotations vs without
#             - with interpolation over the background vs without
#       to plot all would currently be eight different plots which is a bit ridiculous!


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("functions to perform")
    group1.add_argument('--poly', dest='poly', default=False, action='store_true',
                        help="Draw polygons on the supernova remnants instead of simply loading from file (default = False)")
    group1.add_argument('--export', dest='export', default=False, action='store_true',
                        help="Save the resulting polygons, source locations, flux densities, etc to SNR_catalogue.npz (NB: will overwrite any existing measurements) (default = False)")
    group1.add_argument('--plot', dest='makeplots', default=False, action='store_true',
                        help="Make plots of the polygons and sources overlaid on FITS images (NB: will overwrite any existing plots) (default = False)")
    group1.add_argument('--fluxfit', dest='fluxfit', default=False, action='store_true',
                        help="Fit SEDs to the SNRs and make spectral plots (default=False)")
#    group1.add_argument('--fit', dest='fit', default=False, action='store_true',
#                        help="Fit the supernova remnants instead of simply loading from file (default = False)")
#    group1.add_argument("--xm", dest='xm', type=str, default=None,
#                        help='A .fits binary or VO table. The crossmatch between the reference and source catalogue.')
#    group1.add_argument("--infits", dest='infits', type=str, default=None,
#                        help="The fits image(s) to be corrected; enclose in quotes for globbing.")
#    group1.add_argument("--suffix", dest='suffix', type=str, default=None,
#                        help="The suffix to append to rename the output (corrected) fits image(s); e.g., specifying \"warp\" will result in an image like image_warp.fits (no default; if not supplied, no correction will be performed).")
#    group2 = parser.add_argument_group("catalog column names")
#    group2.add_argument("--ra1", dest='ra1', type=str, default='ra',
#                        help="The column name for ra  (degrees) for source catalogue.")
#    group2.add_argument("--dec1", dest='dec1', type=str, default='dec',
#                        help="The column name for dec (degrees) for source catalogue.")
#    group2.add_argument("--ra2", dest='ra2', type=str, default='RAJ2000',
#                        help="The column name for ra  (degrees) for reference catalogue.")
#    group2.add_argument("--dec2", dest='dec2', type=str, default='DEJ2000',
#                        help="The column name for dec (degrees) for reference catalogue.")
#    group3 = parser.add_argument_group("Other")
#    group3.add_argument('--smooth', dest='smooth', default=300.0, type=float,
#                        help="Smoothness parameter to give to the radial basis function (default = 300 pix)")
# Need to update the testing
#    group3.add_argument('--test', dest='test', default=False, action='store_true',
#                        help="run test")

    options = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
 
    # Load the snr info from the files
    snrs = read_snrs()
    if options.poly:
    # Do all the interactive plotting and fitting
        snrs = poly_plots(snrs,options.export)
    else:
        snrs = import_snrs()
    if options.fluxfit:
        snrs = fit_fluxes(snrs,options.export)

    if options.makeplots is True:
        make_plots(snrs)


# Now to add:
# - Make a test simulation fits file that anyone can play with and verify that the code works
# - plotting the final results:
#    - each file as an individual frame with the boundaries and sources shown
#    - if rgb is selected, plot the final RGB image with just the boundary and sources and none of the interpolation boundaries shown
# - saving the output flux densities -- maybe make a new object? class SNR? -- but do I ever want to read in JUST the fluxes?
# - potentially fitting the spectra of sources across the images rather than just a peak subtraction, but that may be too complicated / unnecessary -- need to investigate
# Would also like to add functionality to read off the "z" values for each of the RGB layers I will be plotting, so I can instantly fit the flux of sources and save it somewhere.
