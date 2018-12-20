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

from astropy.visualization import PercentileInterval
#from astropy.visualization import AsinhStretch

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
from defs import normalize

# Normalisation frequency for all functions
normfreq = 200.e6 # 200 MHz

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
    try:
        fit = leastsq(errfunc, pinit, args=(freq_array, flux_array, flux_errors), full_output=1)
    except TypeError:
        print "Could not perform fitting!"
        fit = None
    if fit is not None:
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
       try:
           plt.get_current_fig_manager().window.set_position("GTK_WIN_POS_CENTER")
           plt.get_current_fig_manager().window.resize(1000,500)
       except:
           print "GTK manager not being used, can't put window in the centre."
       return hfig

def poly_plots(snrs):
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

#        rgb = np.dstack([red_data,green_data,blue_data])

        xmax = white[0].header["NAXIS1"]
        ymax = white[0].header["NAXIS2"]

        w = wcs.WCS(white[0].header)
        try:
            pix2deg = white[0].header["CD2_2"]
        except KeyError:
            pix2deg = white[0].header["CDELT2"]
        # Divide by maximum value to normalise
#        rgb = normalize(rgb, np.nanmin(rgb), np.nanmax(rgb))

    # Percentile interval for image normalisation
        pct = 97.0
        interval = PercentileInterval(pct)
#        vmin, vmax = interval.get_limits(data)
#        stretch = AsinhStretch(a=0.1)

        i = interval.get_limits(red_data)
        r = normalize(red_data, *i)
        i = interval.get_limits(green_data)
        g = normalize(green_data, *i)
        i = interval.get_limits(blue_data)
        b = normalize(blue_data, *i)

        rgb = np.dstack([r,g,b])
        rgb[np.isnan(rgb)]=0.0


        def update(val):
            vmin = svmin.val
            vmax = svmax.val
            img_white.set_clim([svmin.val, svmax.val])
            fig.canvas.draw()
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

#        vmin0 = np.min(white_data)
#        vmax0 = np.max(white_data)
        vmin0, vmax0 = interval.get_limits(white_data)
        # include the slider: cribbed from https://stackoverflow.com/questions/5611805/using-matplotlib-slider-widget-to-change-clim-in-image
        svmin = Slider(axmin, "vmin", 2*vmin0, -2*vmin0, valinit=vmin0)
        svmax = Slider(axmax, "vmax", -2*vmax0, 2*vmax0, valinit=vmax0)

        svmin.on_changed(update)
        svmax.on_changed(update)

        fig.suptitle('Left-click = select source, middle-click = source subtract, right-click = select region to exclude from background calculation')

        centerx, centery = w.wcs_world2pix(snr.loc.fk5.ra.value,snr.loc.fk5.dec.value,0)
        print centerx, centery

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
            break
            
        if len(polypick.points.x):
            sources.x, sources.y = w.wcs_pix2world(zip(polypick.points.x,polypick.points.y),0).transpose()
        if len(polypick.exclude.x):
            exclude.x, exclude.y = w.wcs_pix2world(zip(polypick.exclude.x,polypick.exclude.y),0).transpose()
# NB: overwrites data for any old measurements for this live object
        snr.polygon = polygon
        snr.sources = sources
        snr.exclude = exclude

    return snrs

def update_text_file(snr, outtype):
    if outtype == "latex":
        outformat = "{0} & {1:02.0f}:{2:02.0f}:{3:02.0f} & {4:+02.0f}:{5:02.0f}:{6:02.0f} & ${7:3.1f}\\times{8:3.1f}$ & ${9:3.2f}\pm{10:3.2f}$ & ${11:3.2f}\pm{12:3.2f}$ \\\\ \n"
        outvars = [snr.name, \
                snr.loc.fk5.ra.hms.h, snr.loc.fk5.ra.hms.m, snr.loc.fk5.ra.hms.s, \
                snr.loc.fk5.dec.dms.d, abs(snr.loc.fk5.dec.dms.m), abs(snr.loc.fk5.dec.dms.s), \
                60*snr.maj, 60*snr.min, \
                snr.fit_narrow["flux"], snr.fit_narrow["fluxerr"], snr.fit_narrow["alpha"], snr.fit_narrow["alphaerr"]]
        outputfile = "SNR_flux_densities.tex"
        if not os.path.exists(outputfile):
            with open(outputfile, 'w') as output_file:
               output_file.write("Name & RA & Dec & Size & Flux & $\\alpha$ \\\\ \n")
    elif outtype == "csv":
        outformat = "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9}\n"
        outvars = [snr.name, \
                snr.loc.fk5.ra.value, snr.loc.fk5.dec.value, \
                60*snr.maj, 60*snr.min, \
                snr.fit_narrow["flux"], snr.fit_narrow["fluxerr"], snr.fit_narrow["alpha"], snr.fit_narrow["alphaerr"], snr.fit_narrow["chi2red"]]
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

def do_fit(colors,snr):
    freqs = []
    fluxes = []
    bkg = []
    rms = []
    nbeams = []
    for color in colors:
        if color == "white":
            fitsfile = color+"/"+snr.name+".fits"
        else:
            fitsfile = color+"/rpj/"+snr.name+".fits"
        if os.path.exists(fitsfile):
            hdu = fits.open(fitsfile)
            freqs.append(hdu[0].header["FREQ"])
            final_flux, total_flux, bkg_flux, rms_flux, nbeam = find_fluxes(snr.polygon, snr.sources, snr.exclude, fitsfile)
            fluxes.append(final_flux)
            bkg.append(bkg_flux)
            rms.append(rms_flux)
            nbeams.append(nbeam)
        else:
            print "Unable to fit flux for {0} in {1}".format(snr.name, color)
    fluxes = np.array(fluxes)
    freqs = np.array(freqs)
    logfreqs = np.log([f/normfreq for f in freqs])
    logfluxes = np.log(fluxes)
    fluxerrors = [np.sqrt((0.02*s)**2+n*(r**2)) for s,r,n in zip(fluxes,rms,nbeams)]
    logfluxerrors = np.array([ e / f for e,f in zip (fluxerrors, fluxes)])
    alpha, err_alpha, amp, err_amp, chi2red = fit_spectrum(logfreqs[np.where(fluxes>0)],logfluxes[np.where(fluxes>0)],logfluxerrors[np.where(fluxes>0)])

    return snr, dict(zip(freqs, fluxes)), dict(zip(freqs, bkg)), dict(zip(freqs, rms)), dict(zip(freqs, nbeams)), alpha, err_alpha, amp, err_amp, chi2red

def fit_fluxes(snrs):
    for snr in snrs:
# Fit over widebands
        colors = ["white", "red", "green", "blue"]
        snr, snr.flux_wide, snr.bkg_wide, snr.rms_wide, snr.nbeams_wide, alpha, err_alpha, amp, err_amp, chi2red = do_fit(colors,snr)
        if alpha is not None:
            flux200 = np.exp(amp)
            flux200err = err_amp * flux200
            snr.fit_wide = {"flux" : flux200, "fluxerr" : flux200err, "alpha" : alpha, "alphaerr" : err_alpha, "chi2red" : chi2red}
    # Save the plots
            plot_spectrum(snr, wide=True)

# Fit over narrow bands
        print "Fitting spectrum to flux densities from "+snr.name
        colors = ["072-080MHz", "080-088MHz", "088-095MHz", "095-103MHz", "103-111MHz", "111-118MHz", "118-126MHz", "126-134MHz", "139-147MHz", "147-154MHz", "154-162MHz", "162-170MHz", "170-177MHz", "177-185MHz", "185-193MHz", "193-200MHz", "200-208MHz", "208-216MHz", "216-223MHz", "223-231MHz"]
        snr, snr.flux_narrow, snr.bkg_narrow, snr.rms_narrow, snr.nbeams_narrow, alpha, err_alpha, amp, err_amp, chi2red = do_fit(colors,snr)
        if alpha is not None:
            flux200 = np.exp(amp)
            flux200err = err_amp * flux200
            snr.fit_narrow = {"flux" : flux200, "fluxerr" : flux200err, "alpha" : alpha, "alphaerr" : err_alpha, "chi2red" : chi2red}
            print "Reduced chi2 of fit: {0}".format(chi2red)
            plot_spectrum(snr, wide=False)

        if snr.fit_narrow["alpha"] is None and snr.fit_wide["alpha"] is None:
            print "Spectrum could not be fit for {0}.".format(snr.name)
        else:
            for outtype in ["latex", "csv"]:
                update_text_file(snr, outtype)
    return snrs

def plot_spectrum(snr,wide=False):
    if wide is True:
        # frequencies for plotting
        freqs = snr.flux_wide.keys()
        fluxes = snr.flux_wide.values()
        bkg = snr.bkg_wide.values()
        rms = snr.rms_wide.values()
        nbeams = snr.nbeams_wide.values()
    # Nice feature would be to have shading indicating the range of possible fits
    #    fluxerr = snr.fit["fluxerr"]
    #    alphaerr = snr.fit["alphaerr"]
        if snr.fit_wide is not None:
            flux200 = snr.fit_wide["flux"]
            err_flux200 = snr.fit_wide["fluxerr"]
            alpha = snr.fit_wide["alpha"]
            err_alpha = snr.fit_wide["alphaerr"]
        else:
            flux200 = None
            err_flux200 = None
            alpha = None
            err_alpha = None
    else:
        # frequencies for plotting
        freqs = snr.flux_narrow.keys()
        fluxes = snr.flux_narrow.values()
        bkg = snr.bkg_narrow.values()
        rms = snr.rms_narrow.values()
        nbeams = snr.nbeams_narrow.values()
    # Nice feature would be to have shading indicating the range of possible fits
    #    fluxerr = snr.fit["fluxerr"]
    #    alphaerr = snr.fit["alphaerr"]
        if snr.fit_narrow is not None:
            flux200 = snr.fit_narrow["flux"]
            err_flux200 = snr.fit_narrow["fluxerr"]
            alpha = snr.fit_narrow["alpha"]
            err_alpha = snr.fit_narrow["alphaerr"]
        else:
            flux200 = None
            err_flux200 = None
            alpha = None
            err_alpha = None
    #fluxerrors = [np.sqrt((0.02*s)**2+(r**2)) for s,r in zip(fluxes,rms)]
    fluxerrors = [np.sqrt((0.02*s)**2+n*(r**2)) for s,r,n in zip(fluxes,rms,nbeams)]
    plotfreqs = np.linspace(0.9*np.min(freqs),1.1*np.max(freqs),200)

    print "Spectrum fit, making plot for "+snr.name
    if len(freqs) > 4:
        outpng = ("./spectra/"+snr.name+".png")
    else:
        outpng = ("./spectra/"+snr.name+"_wide.png")
    # Plot
    fig = plt.figure(figsize=(10,5))
    # LaTeX-safe
    fig.suptitle(snr.name.replace("_","\,"))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    for ax in ax1, ax2:
        if flux200 is not None:
            ax.plot(plotfreqs, powerlaw(plotfreqs, flux200/(normfreq**alpha), alpha),label="$\\alpha={0:3.2f}\pm{1:3.2f}$ ".format(alpha, err_alpha))     # Scipy Fit
        ax.legend()
        ax.errorbar(freqs, bkg, yerr=rms, fmt = "r.", alpha = 0.5)  # Background
        ax.errorbar(freqs, fluxes, yerr=fluxerrors, fmt = "k.") # Data
        ax.set_ylabel("S / Jy")
        ax.set_xlabel("Frequency / MHz")
    ax2.set_xlim(0.8*np.min(freqs), 1.2*np.max(freqs))
    ax2.set_ylim(0.3*np.min([np.min(fluxes),np.min(bkg)]),3*np.max([np.max(fluxes),np.max(bkg)]))
    if os.path.exists(outpng):
        renumber(outpng)
    fig.savefig(outpng)
    fig.clf()

# Naming convention for outputfiles
# types/<snr.name><.type><.ii>
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

def filter_latitude(snrs):
    filtered = []
    for snr in snrs:
        l = snr.loc.galactic.l.value
        if (((l>180) and (l<240)) or (l>340) or (l<60)):
            filtered.append(snr)
    return filtered
       
def test_export_snr(snrs):
    tempsnrs = []
    for snr in snrs:
        file_ext = "pkl"
# Temporary HACK
        output_dir = "old_pkls/"
        output_file = "{0}/{1}.{2}".format(output_dir, snr.name, file_ext)
        if os.path.exists(output_file):
            print "{0} has already been written.".format(snr.name)
        else:
            print "{0} will be measured.".format(snr.name)
            tempsnrs.append(snr)
    return tempsnrs

def export_snrs(snrs):
    for snr in snrs:
        file_ext = "pkl"
        output_dir = "pkls/"
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
        if snr.fit_narrow is None:
            if snr.fit_wide is None:
                fittype = "none"
                sf = ""
                err_sf = ""
                sa = ""
                err_sa = ""
                sc = ""
            else:
                fittype = "wide"
                sf = snr.fit_wide["flux"]
                err_sf = snr.fit_wide["fluxerr"]
                sa = snr.fit_wide["alpha"]
                err_sa = snr.fit_wide["alphaerr"]
                sc = snr.fit_wide["chi2red"]
        else:
            fittype = "narrow"
            sf = snr.fit_narrow["flux"]
            err_sf = snr.fit_narrow["fluxerr"]
            sa = snr.fit_narrow["alpha"]
            err_sa = snr.fit_narrow["alphaerr"]
            sc = snr.fit_narrow["chi2red"]
        t = Table([[snr.name], [snr.loc.fk5.ra.value], [snr.loc.fk5.dec.value], [snr.maj], [snr.min], [snr.pa], [snr.known], [sf], [err_sf], [sa], [err_sa], [sc], [fittype]], names=("Name", "RAJ2000", "DEJ2000", "major", "minor", "pa", "known", "flux_fitted", "err_flux_fitted", "alpha", "err_alpha", "chi2red", "fittype"))
        if snr.flux_narrow is not None:
            for key, value in snr.flux_narrow.iteritems():
                col = Column(name = "flux_{0}".format(key), data = [value])
                t.add_column(col)
        if snr.flux_wide is not None:
            for key, value in snr.flux_wide.iteritems():
                col = Column(name = "flux_{0}".format(key), data = [value])
                t.add_column(col)
        
        t.write(output_file, format=file_ext)

def import_old_snr(snr):
    print "Checking for existing measurements for {0}".format(snr.name)
    file_ext = "pkl"
#HACK
    input_dir = "/home/tash/Dropbox/MWA/SNR_search/old_pkls"
    input_file = "{0}/{1}.{2}".format(input_dir, snr.name, file_ext)
#HACK to preserve RA and Dec
    coords = snr.loc
    if os.path.exists(input_file):
        ifile = open(input_file, "rb")
        snr = pickle.load(ifile)
        snr.loc = coords
    else:
        print "No existing measurements for {0}".format(snr.name)

    new_snr = SNR()
    new_snr.name = snr.name
    new_snr.loc = snr.loc
    new_snr.maj = snr.maj
    new_snr.min = snr.min
    new_snr.known = snr.known
    new_snr.polygon = snr.polygon
    new_snr.sources = snr.sources
    new_snr.exclude = snr.exclude
    if len(snr.flux.keys()) > 4:
        new_snr.flux_narrow = snr.flux
        new_snr.bkg_narrow = snr.bkg
        new_snr.rms_narrow = snr.rms
        new_snr.nbeams_narrow = snr.nbeams
        new_snr.fit_narrow = snr.fit
# NB: hopefully I won't need this because I will regenerate everything right away
#        new_snr.fit_narrow["flux"] = snr.fit["flux"] * (200./150.)**snr.fit["alpha"]
#        new_snr.fit_narrow["fluxerr"] = snr.fit["fluxerr"] * (new_snr.fit_narrow["flux"]/snr.fit["flux"])
#        new_snr.fit_narrow["alpha"] = snr.fit["alpha"]
#        new_snr.fit_narrow["err_alpha"] = snr.fit["err_alpha"]
#        new_snr.fit_narrow["chi2red"] = snr.fit["chi2red"]
    else:
        new_snr.flux_wide = snr.flux
        new_snr.bkg_wide = snr.bkg
        new_snr.rms_wide = snr.rms
        new_snr.nbeams_wide = snr.nbeams
        new_snr.fit_wide = snr.fit
    return new_snr

def import_snr(snr):
    print "Checking for existing measurements for {0}".format(snr.name)
    file_ext = "pkl"
#HACK to direct path -- would be nice to tidy this
    input_dir = "/home/tash/Dropbox/MWA/SNR_search/pkls"
    input_file = "{0}/{1}.{2}".format(input_dir, snr.name, file_ext)
# Preserve Name, RA, Dec, a, b, pa from candidates.txt
    coords = snr.loc
    snrmaj = snr.maj
    snrmin = snr.min
    snrpa = snr.pa
    snrname = snr.name
    if os.path.exists(input_file):
        ifile = open(input_file, "rb")
        snr = pickle.load(ifile)
        snr.loc = coords
        snr.maj = snrmaj
        snr.min = snrmin
        snr.pa = snrpa
        snr.name = snrname
#    else:
#        if abs(snr.loc.galactic.b)<0.001:
# Check if it's a +/- 0.0 SNR
#        input_file = "{0}/{1}{2}.{3}".format(input_dir, snr.loc.galactic.l, snr.loc.galactic.b, file_ext)
    else:
        print "No existing measurements for {0}".format(snr.name)
    return snr

def make_plots(snrs):
    for snr in snrs:
        print "Making attractive FITS image plot for "+snr.name
# Load image data
        mask = fits.open("white/"+snr.name+"_mask.fits")
        white = fits.open("white/"+snr.name+".fits")
        red = fits.open("red/rpj/"+snr.name+".fits")
        green = fits.open("green/rpj/"+snr.name+".fits")
        blue = fits.open("blue/rpj/"+snr.name+".fits")
        mask_data = mask[0].data
# Set mask data to 1.0, NaNs will not be plotted
        mask_data[np.where(np.logical_not(np.isnan(mask_data)))] = 1.0
        white_data = white[0].data
        red_data = red[0].data
        green_data = green[0].data
        blue_data = blue[0].data
        rgb = np.dstack([red_data,green_data,blue_data])
        rgb = normalize(rgb, np.nanmin(rgb), np.nanmax(rgb))
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

# Using http://docs.astropy.org/en/stable/visualization/wcsaxes/index.html
        fig = plt.figure(figsize=(6,6))
        fig.set_tight_layout(True)
        ax_white = fig.add_subplot(221, projection=w)
        ax_white_mask = fig.add_subplot(222, projection=w)
        ax_rgb = fig.add_subplot(223, projection=w)
        ax_rgb_mask = fig.add_subplot(224, projection=w)
        img_white = ax_white.imshow(white_data, cmap="gray",origin="lower")
        img_white_mask = ax_white_mask.imshow(white_data, cmap="gray",origin="lower")
        img_white_mask = ax_white_mask.imshow(mask_data, cmap="Blues",origin="lower", alpha = 0.2)
        img_rgb = ax_rgb.imshow(rgb,origin="lower")
        img_rgb_mask = ax_rgb_mask.imshow(rgb,origin="lower")
        img_rgb_mask = ax_rgb_mask.imshow(mask_data, cmap="Blues",origin="lower", alpha = 0.2)
        for ax in ax_white_mask, ax_rgb_mask:
            if len(local_polygon.x):
                ax.plot(local_polygon.x,local_polygon.y,**restsnr)
            if len(local_sources.x):
                ax.plot(local_sources.x,local_sources.y,**srcmark)
            if len(local_exclude.x):
                ax.plot(local_exclude.x,local_exclude.y,**restexc)
# Tedious removal of axis labels
        axes = [ax_white, ax_white_mask, ax_rgb, ax_rgb_mask]
        latleft = [ True, False, True, False ]
        latright = [ False, True, False, True ]
        lonbottom = [ False, False, True, True ]
        lontop = [ True, True, False, False ]
        for ax, lal, lar, lob, lot in zip(axes, latleft, latright, lonbottom, lontop):
            ax.set_xlim(0,xmax)
            ax.set_ylim(0,ymax)
            lon = ax.coords['ra']
            lat = ax.coords['dec']
            if lob:
                lon.set_axislabel("Right Ascension (J2000)")
                lon.set_major_formatter('hh:mm')
            if lal:
                lat.set_axislabel("Declination (J2000)")
                lat.set_major_formatter('dd:mm')
            lon.set_ticklabel_visible(lob)
            lat.set_ticklabel_visible(lal)
            lon.set_ticks_visible(lob)
            lat.set_ticks_visible(lal)
            overlay = ax.get_coords_overlay("fk5")
            overlay.grid(axes = ax, color='white', ls='dotted')
            overlay["ra"].set_ticklabel_visible(lot)
            if lot:
                overlay["ra"].set_major_formatter('hh:mm')
            overlay["dec"].set_ticklabel_visible(lar)
            if lar:
                overlay["dec"].set_major_formatter('hh:mm')

        output_file = "plots/"+snr.name+".png"
        if os.path.exists(output_file):
            renumber(output_file)
        fig.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)
#        fig2.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
#        fig1.savefig("int_flux_ratio.png",pad_inches=0.1,bbox_inches="tight")
#        fig2.savefig("divergence.png",pad_inches=0.1,bbox_inches="tight")

        fig.savefig(output_file,pad_inches=0.1,bbox_inches="tight")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("functions to perform")
    group1.add_argument('--poly', dest='poly', default=False, action='store_true',
                        help="Draw polygons on the supernova remnants instead of simply loading from file (default = False)")
    group1.add_argument('--overdraw', dest='overdraw', default=False, action='store_true',
                        help="Draw polygons even if a SNR has already been measured (default = False)")
    group1.add_argument('--fluxfit', dest='fluxfit', default=False, action='store_true',
                        help="Subtract backgrounds and sources, fit SEDs to the SNR flux densities, and make spectral plots (default=False)")
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

# Filter those outside of my search range

    snrs = filter_latitude(snrs)

    updated_snrs = []
    for snr in snrs:
# HACK to import old style snr with 150-MHz centers and narrow OR wide flux density measurements
#        snr = import_old_snr(snr)
        snr = import_snr(snr)
        updated_snrs.append(snr)
    snrs = updated_snrs

    print "Read {0} SNRs".format(len(snrs))
    # Do all the interactive plotting and fitting
    if options.poly:
    # Don't overwrite polygons for already-measured SNRs
        if not options.overdraw:
            print "Overdraw not selected: filtering out SNRs which have already been fitted."
            print "Rerun without --poly if you wanted to re-measure their flux densities without drawing."
            snrs = test_export_snr(snrs)
    # Only snrs which have been fitted will be returned by the function
        print "Measuring {0} SNRs".format(len(snrs))
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
