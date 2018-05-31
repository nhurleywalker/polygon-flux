#!/usr/bin/env python

__author__ = "Natasha Hurley-Walker"
__date__ = "21/05/2018"

import os
import sys
import shutil
import glob

import numpy as np
import pickle
from astropy.io import fits
from astropy import wcs
from astropy.table import Table, Column

# Need a least-squares estimator that gives a useable error estimate
from scipy.optimize import leastsq

import matplotlib
#matplotlib.use('Agg') # So does not use display -- only good if just making plots
from matplotlib.widgets import Slider, RadioButtons
from matplotlib import pyplot as plt

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

def poly_plots(snrs):
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
# NB: overwrites data for any old measurements for this live object
        snr.polygon = polygon
        snr.sources = sources
        snr.exclude = exclude
        snrs2.append(snr)

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

def do_fit(colors,snr,normfreq):
    freqs = []
    fluxes = []
    bkg = []
    for color in colors:
        if color == "white":
            fitsfile = color+"/"+snr.name+".fits"
        else:
            fitsfile = color+"/rpj/"+snr.name+".fits"
        hdu = fits.open(fitsfile)
        freqs.append(hdu[0].header["FREQ"])
        final_flux, total_flux, source_flux, bkg_flux = find_fluxes(snr.polygon, snr.sources, snr.exclude, fitsfile)
        fluxes.append(final_flux)
        bkg.append(bkg_flux)

    logfreqs = np.log([f/normfreq for f in freqs])
    logfluxes = np.log(fluxes)

    print freqs
    print fluxes
# 10% errors to start with
    fluxerrors = [0.1*s for s in fluxes]
    alpha, err_alpha, amp, err_amp, chi2red = fit_spectrum(logfreqs,logfluxes,0.1)
# Save to SNR object for later
    snr.flux = dict(zip(freqs, fluxes))
    snr.bkg = dict(zip(freqs, bkg))

    return snr, alpha, err_alpha, amp, err_amp, chi2red

def fit_fluxes(snrs):
    normfreq = 150000000. # 150 MHz
    for snr in snrs:
        print "Fitting spectrum to flux densities from "+snr.name
        colors = ["072-080MHz", "080-088MHz", "088-095MHz", "095-103MHz", "103-111MHz", "111-118MHz", "118-126MHz", "126-134MHz", "139-147MHz", "147-154MHz", "154-162MHz", "162-170MHz", "170-177MHz", "177-185MHz", "185-193MHz", "193-200MHz", "200-208MHz", "208-216MHz", "216-223MHz", "223-231MHz"]
        snr, alpha, err_alpha, amp, err_amp, chi2red = do_fit(colors,snr,normfreq)

        print "Reduced chi2 of fit: {0}".format(chi2red)
# If fitting across all the sub-bands failed, use just the wideband images
#        if not alpha or chi2red > 1.93:
#            colors = ["white", "red", "green", "blue"]
#            print "Sub-band fit failed: using wideband images."
#            snr, alpha, err_alpha, amp, err_amp, chi2red = do_fit(colors,snr,normfreq)

        # frequencies for plotting
        freqs = snr.flux.keys()
        fluxes = snr.flux.values()
        bkg = snr.bkg.values()
        fluxerrors = [0.1*s for s in fluxes]
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
            ax1.scatter(freqs, bkg, marker="o", c = "red")  # Background
            ax1.set_ylabel("S / Jy")
            ax1.set_xlabel("Frequency / MHz")
            ax1.legend()
            ax2=example.add_subplot(1,2,2)
            ax2.set_xscale("log")
            ax2.set_yscale("log")
            ax2.plot(plotfreqs, powerlaw(plotfreqs, flux150/(normfreq**alpha), alpha),label="alpha=${0:3.1f}$ ".format(alpha))     # Scipy Fit
            ax2.errorbar(freqs, fluxes, yerr=fluxerrors, fmt='k.')  # Data
            ax2.scatter(freqs, bkg, marker="o", c = "red")  # Background
    #        ax2.set_xlim(left=xmin01,right=xmax01)
    #       ax2.set_ylim([min(np.exp(flux_array)),max(np.exp(flux_array))])
            ax2.set_ylabel("S / Jy")
            ax2.set_xlabel("Frequency / MHz")
            ax2.legend()
            if os.path.exists(outpng):
                renumber(outpng)
            example.savefig(outpng)
            example.clf()
        
            snr.fit = {"flux" : flux150, "fluxerr" : flux150err, "alpha" : alpha, "alphaerr" : err_alpha, "chi2red" : chi2red}
# Add table output
            for outtype in ["latex", "csv"]:
                update_text_file(snr, outtype)
    return snrs

# Naming convention for outputfiles
# npzs/<snr.name><.npz><.ii>
# where ii is a 2-digit integer
def renumber(output_file):
   # name, file_ext = os.path.splitext(output_file)
    allfiles = sorted(glob.glob(output_file+"*"), reverse=True)
    for fl in allfiles:
       iis = fl.split(".")[-1]
       try:
           ii = int(iis)
           ii += 1
           pn = os.path.splitext(fl)[0]
       except:
           ii = 0
           pn = fl
       shutil.move(fl, "{0}.{1:02d}".format(pn,ii))

def export_snrs(snrs):
    print "Exporting SNR catalogue"
    file_ext = "pkl"
    output_dir = "pkls/"
    for snr in snrs:
        print "Pickling {0}".format(snr.name)
        output_file = "{0}/{1}.{2}".format(output_dir, snr.name, file_ext)
        if os.path.exists(output_file):
            renumber(output_file)
        ofile = open(output_file, "wb")
        pickle.dump(snr, ofile)
        file_ext = "fits"
        output_dir = "cats/"
        output_file = "{0}/{1}.{2}".format(output_dir, snr.name, file_ext)
        print "Writing FITS catalogue for {0}".format(snr.name)
        if os.path.exists(output_file):
            renumber(output_file)
# FINISH THIS
        if snr.fit is None:
            sf = ""
            sa = ""
            sc = ""
        else:
            sf = snr.fit["flux"]
            sa = snr.fit["alpha"]
            sc = snr.fit["chi2red"]
        t = Table([[snr.name], [snr.loc.ra.value], [snr.loc.dec.value], [snr.maj], [snr.min], [snr.known], [sf], [sa], [sc]], names=("Name", "RAJ2000", "DEJ2000", "major", "minor", "known", "flux_fitted", "alpha", "chi2red"))
        print t
        print snr.flux
        if snr.flux is not None:
            for key, value in snr.flux.iteritems():
                print key, value
                col = Column(name = "flux_{0}".format(key), data = [value])
                print col
                t.add_column(col)
    # Still to output, probably in a different table: self.flux = None # A dictionary of frequencies to flux densities
        t.write(output_file, format=file_ext)

def import_snr(snr):
    print "Checking for existing measurements for {0}".format(snr.name)
    file_ext = "pkl"
    input_dir = "pkls"
    input_file = "{0}/{1}.{2}".format(input_dir, snr.name, file_ext)
    if os.path.exists(input_file):
        ifile = open(input_file, "rb")
        snr = pickle.load(ifile)
    else:
        print "No existing measurements for {0}".format(snr.name)
    return snr

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
        output_file = "plots/"+snr.name+".png"
        if os.path.exists(output_file):
            renumber(output_file)
        fig.savefig(output_file)

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
    group1.add_argument('--fluxfit', dest='fluxfit', default=False, action='store_true',
                        help="Fit SEDs to the SNRs, create interpolated background images and make spectral plots (default=False)")
    group1.add_argument('--export', dest='export', default=False, action='store_true',
                        help="Save the resulting polygons, source locations, flux densities, etc to npz files (default=False)")
    group1.add_argument('--plot', dest='makeplots', default=False, action='store_true',
                        help="Make plots of the polygons and sources overlaid on FITS images (default = False)")
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
 
    snrs = read_snrs()
    updated_snrs = []
    for snr in snrs:
        snr = import_snr(snr)
        updated_snrs.append(snr)
    snrs = updated_snrs

    # Do all the interactive plotting and fitting
    if options.poly:
    # Only snrs which have been fitted will be returned by the function
        snrs = poly_plots(snrs)

    # Fit flux densities across the band using the measured polygons
    if options.fluxfit:
        fitsnrs = []
        # Necessary in case you load these from existing files which somehow don't have polygons
        for snr in snrs:
            if snr.polygon is not None:
                fitsnrs.append(snr)
            else:
                print "No measured polygons for {0}, cannot fit fluxes".format(snr.name)
        snrs = fit_fluxes(fitsnrs)

    if options.export:
        export_snrs(snrs)

    if options.makeplots is True:
        plotsnrs = []
        for snr in snrs:
            if snr.polygon is not None:
                plotsnrs.append(snr)
            else:
                print "No measured polygons for {0}, cannot plot".format(snr.name)
        make_plots(plotsnrs)

# Now to add:
# - Make a test simulation fits file that anyone can play with and verify that the code works
# - plotting the final results:
#    - each file as an individual frame with the boundaries and sources shown
#    - if rgb is selected, plot the final RGB image with just the boundary and sources and none of the interpolation boundaries shown
# - saving the output flux densities -- maybe make a new object? class SNR? -- but do I ever want to read in JUST the fluxes?
# - potentially fitting the spectra of sources across the images rather than just a peak subtraction, but that may be too complicated / unnecessary -- need to investigate
# Would also like to add functionality to read off the "z" values for each of the RGB layers I will be plotting, so I can instantly fit the flux of sources and save it somewhere.
