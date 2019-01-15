#!/usr/bin/env python

__author__ = "Natasha Hurley-Walker"
__date__ = "09/01/2018"

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
from astropy import units as u
from astropy.nddata import Cutout2D

from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve

from astropy.io.votable import parse_single_table
from astropy.visualization import PercentileInterval

from astropy.table import Column
from astropy.table import Table

from decimal import Decimal

import numpy as np

import defs
from defs import SNR
from defs import read_snrs
from measure_snr import import_snr
from measure_snr import fit_spectrum
from measure_snr import powerlaw
import integrated_polys

import table_formatting

normfreq = 200.e6

#def combine_cols(data, errors, name):
## First loop identifies the largest number of decimal places I need to show
#    mine = 0
#    for i in errors:
#        if i is not None:
#            e = Decimal("{0:2.1g}".format(i)).as_tuple().exponent
#            if e < mine:
#                mine = e
#    if mine < 0:
#        fmt = "g"
##    else:
##   Need to figure out how to format in the case of a load of errors > 0, to make this function general-purpose

## Second loop does the formating
#    rs = []
#    for i, j in zip(data, errors):
#        if i is not None and j is not None:
#            r = ("${0:2.{1}f}\pm{2:2.{1}f}$".format(i, abs(mine), j))
#        else:
#            r = ("--")
#        rs.append(r)
#
#    newcol = Column(rs, name=name)
#    return newcol

def normalize(arr, vmin, vmax):
    nor = (arr - vmin) / (vmax - vmin)
    nor[np.where(nor<0.0)] = 0.0
    nor[np.where(nor>1.0)] = 1.0
    return nor

def regrid(fitsimage, outname, template, overwrite=False):
    if not os.path.exists(outname) or overwrite is True:
         os.system("fits op=xyin in="+fitsimage+" out=gleam.im")
         os.system("fits op=xyin in="+template+" out=temp.im")
         os.system("regrid in=gleam.im out=gleam_rg.im tin=temp.im tol=0.0 axes=1,2")
         os.system("fits op=xyout in=gleam_rg.im out="+outname)
         shutil.rmtree("gleam.im")
         shutil.rmtree("temp.im")
         shutil.rmtree("gleam_rg.im")

def cand_out(ra, dec, diam):
# Export a fake candidates.txt so I don't have to rewrite my SNR reading code
    outputfile = "candidates.txt"
    if os.path.exists(outputfile):
        os.remove(outputfile)
    for r, d, di in zip(ra, dec, diam):
        outformat = "{0} {1} {2}\n"
        outvars = [r, d, di/60.]
        with open(outputfile, 'a') as output_file:
           output_file.write(outformat.format(*outvars))

def find_sources(snr=None):
    catalogue = "/media/data/MWA/GLEAM/GP/GP_Catalogue_Paper/GLEAM_GAL_joined.fits"

    hdu_cat = fits.open(catalogue)
    cat = hdu_cat[1].data

    catalog = SkyCoord(cat["RAJ2000"]*u.deg, cat["DEJ2000"]*u.deg)
    c = SkyCoord([snr.loc.ra], [snr.loc.dec], frame="fk5", unit = (u.deg, u.deg))
    idx, idxcatalog, d2d, d3d = catalog.search_around_sky(c, (snr.maj/2.)*u.deg)

    if len(idxcatalog > 0):
        sources = cat[idxcatalog]
# Filter out all sources with poorly-fit spectral indices
        print "total number of sources:"
        print len(sources)
        sources = sources[np.where(np.logical_not(np.isnan(sources["alpha"])))]
        print "sources with well-fit spectral indices:"
        print len(sources)
    else:
        sources = None

    # Integrate flux density of sources for print statement only
    total_flux = 0.0
    alphas = []
    for source in sources:
        total_flux += source["int_flux_fit_200"]
        print source["Name"], source["int_flux_fit_200"], source["alpha"], total_flux
        alphas.append(source["alpha"])
    print "Subtracting $S_\mathrm{{200MHz}}={0}$\,Jy, with median $\\alpha={1}$".format(total_flux, np.median(alphas))

    return sources


def refit(snr=None, high=False, sourcesub=False):
    freqs = snr.flux_narrow.keys()
    fluxes = snr.flux_narrow.values()
    bkg = snr.bkg_narrow.values()
    rms = snr.rms_narrow.values()
    nbeams = snr.nbeams_narrow.values()
    fluxerrors = [np.sqrt((0.02*s)**2+n*(r**2)) for s,r,n in zip(fluxes,rms,nbeams)]

# Only want to select frequencies above 150MHz
    if high is True:

        high_freqs = []
        high_fluxes = []
        high_bkg = []
        high_rms = []
        high_nbeams = []
        high_fluxerrors = []

        for f, s, b, r, n, e in zip(freqs, fluxes, bkg, rms, nbeams, fluxerrors):
            if f > 150.e6:
                high_freqs.append(f)
                high_fluxes.append(s)
                high_bkg.append(b)
                high_rms.append(r)
                high_nbeams.append(n)
                high_fluxerrors.append(e)
    else:
        high_freqs = freqs
        high_fluxes = fluxes
        high_bkg = bkg
        high_rms = rms
        high_nbeams = nbeams
        high_fluxerrors = fluxerrors

    if sourcesub is True:
        sources = find_sources(snr)
        if sources is not None:
            src_fluxes = []
            for freq in high_freqs:
                src_flux = 0.0
                for source in sources:
                    src_flux += source["int_flux_fit_200"]*(freq/normfreq)**source["alpha"]
                src_fluxes.append(src_flux)
            high_fluxes = np.array(high_fluxes) - np.array(src_fluxes)
        else:
            print "No sources detected for "+snr.name
            src_fluxes = None

    bkgfreqs = list(freqs)

    high_fluxes = np.array(high_fluxes)
    high_freqs = np.array(high_freqs)
    high_logfreqs = np.log([f/normfreq for f in high_freqs])
    high_logfluxes = np.log(high_fluxes)
    high_logfluxerrors = np.array([ e / f for e,f in zip (high_fluxerrors, high_fluxes)])

    print high_fluxes, high_freqs

    plotfreqs = np.linspace(0.9*np.min(freqs),1.1*np.max(freqs),200)

    alpha, err_alpha, amp, err_amp, chi2red = fit_spectrum(high_logfreqs[np.where(high_fluxes>0)],high_logfluxes[np.where(high_fluxes>0)],high_logfluxerrors[np.where(high_fluxes>0)])
    flux200 = np.exp(amp)
    flux200err = err_amp * flux200

    print "Spectrum fit, making plot for "+snr.name
    outpng = (snr.name+"_spectrum.png")
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
            ax.plot(plotfreqs, powerlaw(plotfreqs, flux200/(normfreq**alpha), alpha),label="alpha=${0:3.2f}$ ".format(alpha))     # Scipy Fit
        ax.legend()
        ax.errorbar(bkgfreqs, bkg, yerr=rms, fmt = "r.", alpha = 0.5)  # Background
        ax.errorbar(freqs, fluxes, yerr=fluxerrors, fmt = "k.", alpha = 0.3) # All Data
        ax.errorbar(high_freqs, high_fluxes, yerr=high_fluxerrors, fmt = "k.") # Fitted Data
        if sourcesub is True:
            ax.scatter(high_freqs, np.array(high_fluxes)+np.array(src_fluxes), color="yellow", marker=".", alpha = 0.5)
            ax.scatter(high_freqs, src_fluxes, color="green", marker=".", alpha = 0.5)
        ax.set_ylabel("S / Jy")
        ax.set_xlabel("Frequency / MHz")
    ax2.set_xlim(0.8*np.min(freqs), 1.2*np.max(freqs))
    ax2.set_ylim(0.3*np.min([np.min(fluxes),np.min(bkg)]),3*np.max([np.max(fluxes),np.max(bkg)]))
    fig.savefig(outpng)
    fig.clf()

    return alpha, err_alpha, flux200, flux200err

# read from candidates.txt
snrs = read_snrs()
# Get measurements from the PKL files
updated_snrs = []
for snr in snrs:
    snr = import_snr(snr)
    updated_snrs.append(snr)
snrs = updated_snrs

## How to fit S_200 and alpha for SNRs

# SNRs where we first need to subtract the GLEAM sources within the shell
sub_snrs = ["SNR_G2.2+2.8", "SNR_G25.4-1.9", "SNR_G38.7-1.1", "SNR_G230.5+1.3"]

# SNRs for which we have successfully fit the full spectrum
good_snrs = ["SNR_G355.4+2.8", "SNR_G356.5-1.9", "SNR_G18.9-1.3", "SNR_G28.4+0.2", "SNR_G0.2-9.7", "SNR_G35.4-0.0", "SNR_G230.5+1.3"]

# SNRs for which we can fit using a combination of our data and Effelsberg
E11_snrs = ["SNR_G24.1-0.3", "SNR_G232.2+2.1", "SNR_G19.2-3.1", "SNR_G2.2+2.8", "SNR_G25.3-1.9", "SNR_G38.7-1.1", "SNR_G21.8+0.2", "SNR_G23.1+0.2", "SNR_G25.4-1.9", "SNR_G19.7-0.7", "SNR_G28.8-0.5"]

# SNRs for which we can fit using a combination of our data and Molonglo GPS
MGPS_snrs = ["SNR_G351.9+0.2", "SNR_G349.1-0.8", "SNR_G351.5+0.2", "SNR_G351.0-0.6"]

# SNRs which we can only fit using data points > 150MHz due to non-thermal contamination
high_snrs = ["SNR_G350.8+0.7", "SNR_G353.1+0.8", "SNR_G358.4-0.8"]

## Classes of SNRs

c1_snrs = ["SNR_G232.2+2.1", "SNR_G355.4+2.8", "SNR_G351.9+0.2", "SNR_G356.5-1.9", "SNR_G19.7-0.7", "SNR_G18.9-1.3", "SNR_G21.8+0.2", "SNR_G23.1+0.2", "SNR_G24.1-0.3", "SNR_G19.2-3.1", "SNR_G28.4+0.2", "SNR_G28.8-0.5", "SNR_G25.4-1.9"]

c2_snrs = ["SNR_G230.5+1.3", "SNR_G350.8+0.7", "SNR_G349.1-0.8", "SNR_G351.5+0.2", "SNR_G353.1+0.8", "SNR_G351.0-0.6", "SNR_G2.2+2.8", "SNR_G35.4-0.0", "SNR_G38.7-1.1"]

c3_snrs = ["SNR_G358.4-0.8", "SNR_G0.2-9.7", "SNR_G20.2-0.3"]

## Morphology of SNRs

full_snrs = ["SNR_G232.2+2.1", "SNR_G355.4+2.8", "SNR_G351.9+0.2", "SNR_G356.5-1.9", "SNR_G19.7-0.7", "SNR_G18.9-1.3", "SNR_G21.8+0.2", "SNR_G23.1+0.2", "SNR_G24.1-0.3", "SNR_G19.2-3.1", "SNR_G28.4+0.2", "SNR_G28.8-0.5", "SNR_G230.5+1.3", "SNR_G349.1-0.8", "SNR_G2.2+2.8", "SNR_G25.4-1.9", "SNR_G38.7-1.1"]

part_snrs = ["SNR_G350.8+0.7", "SNR_G351.5+0.2", "SNR_G353.1+0.8", "SNR_G351.0-0.6", "SNR_G35.4-0.0", "SNR_G358.4-0.8", "SNR_G0.2-9.7", "SNR_G20.2-0.3"]


# Create master "properties" VO table, that will later be converted to LaTeX
# columns = Name, RA, Dec, a, b, PA, S200, alpha
# Morphology, Class, and Notes, are defined by hand in Overleaf

Names = []
ls = []
bs = []
RAs = []
Decs = []
majs = []
mins = []
pas = []
S200s = []
err_S200s = []
alphas = []
err_alphas = []
classes = []
morphs = []
ancs = []

for snr in snrs:
    print snr.name
    Names.append(snr.name.replace("SNR_","").replace("G","G\,"))
    ls.append(snr.loc.galactic.l.value)
    bs.append(snr.loc.galactic.b.value)
# Change to sexagecimal
    RAs.append(snr.loc.fk5.ra.to_string(unit=u.hour, sep=(" ", " ", " "), pad=True, precision=0))
    Decs.append(snr.loc.fk5.dec.to_string(unit=u.deg, sep=(" ", " ", " "), pad=True, precision=0, fields =2))
    majs.append("{0:3.0f}".format(60*snr.maj))
    mins.append("{0:3.0f}".format(60*snr.min))
    pas.append("{0:3.0f}".format(snr.pa))
    sourcesub = False
    if snr.name in sub_snrs:
        sourcesub = True
    if snr.name in good_snrs:
        ancs.append("--")
        if sourcesub is False:
            S200s.append(snr.fit_narrow["flux"])
            err_S200s.append(snr.fit_narrow["fluxerr"])
            alphas.append(snr.fit_narrow["alpha"])
            err_alphas.append(snr.fit_narrow["alphaerr"])
        else:
            alpha, err_alpha, flux200, flux200err = refit(snr=snr, high=False, sourcesub=sourcesub)
            S200s.append(flux200)
            err_S200s.append(flux200err)
            alphas.append(alpha)
            err_alphas.append(err_alpha)
    elif snr.name in high_snrs:
# Do the fitting to get these numbers
        alpha, err_alpha, flux200, flux200err = refit(snr=snr, high=True, sourcesub=sourcesub)
        S200s.append(flux200)
        err_S200s.append(flux200err)
        alphas.append(alpha)
        err_alphas.append(err_alpha)
        ancs.append("--")
    elif snr.name in E11_snrs:
# Get the value from Effelsberg
        S_11, tot_11, bkg_11, rms_11, nbeam_11 = integrated_polys.measure_E11(snr)
        snr.flux_narrow[2695000000] = S_11
        snr.bkg_narrow[2695000000] = bkg_11
        snr.rms_narrow[2695000000] = rms_11
        snr.nbeams_narrow[2695000000] = nbeam_11
# Do the fitting to get these numbers
        alpha, err_alpha, flux200, flux200err = refit(snr, high=False, sourcesub=sourcesub)
        S200s.append(flux200)
        err_S200s.append(flux200err)
        alphas.append(alpha)
        err_alphas.append(err_alpha)
        ancs.append("E11")
    elif snr.name in MGPS_snrs:
# Get the value from MGPS
        S_843, tot_843, bkg_843, rms_843, nbeam_843 = integrated_polys.measure_MGPS(snr)
        snr.flux_narrow[843000000] = S_843
        snr.bkg_narrow[843000000] = bkg_843
        snr.rms_narrow[843000000] = rms_843
        snr.nbeams_narrow[843000000] = nbeam_843
# Do the fitting to get these numbers
        alpha, err_alpha, flux200, flux200err = refit(snr, high=False, sourcesub=sourcesub)
        S200s.append(flux200)
        err_S200s.append(flux200err)
        alphas.append(alpha)
        err_alphas.append(err_alpha)
        ancs.append("MGPS")
    else:
        S200s.append(None)
        err_S200s.append(None)
        alphas.append(None)
        err_alphas.append(None)
        ancs.append("--")
# Make class column list
    if snr.name in c1_snrs:
        classes.append("\\textsc{I}")
    elif snr.name in c2_snrs:
        classes.append("\\textsc{II}")
    elif snr.name in c3_snrs:
        classes.append("\\textsc{III}")
# Make morphology column list
    if snr.name in full_snrs:
        morphs.append("Shell")
    elif snr.name in part_snrs:
        morphs.append("Partial shell")


for col in Names, RAs, Decs, majs, mins, pas, S200s, err_S200s, alphas, err_alphas, classes:
    print len(col)

t = Table([Names, ls, bs, RAs, Decs, majs, mins, pas, S200s, err_S200s, alphas, err_alphas], names = ("Name", "GLON", "GLAT", "RA", "Dec", "a", "b", "pa", "S200", "err_S200", "alpha", "err_alpha"))

# Create master "raw measurements" VO table, that will later be converted to LaTeX
# Not yet implemented because it would be 20 x 26 and maybe just not that interesting!!

# Make nice LaTeX-formatted measurement columns
t.add_column(Column(table_formatting.combine_cols(S200s, err_S200s, uniform=False), name="$S_\mathrm{200\,MHz}$"))
t.add_column(Column(table_formatting.combine_cols(alphas, err_alphas, uniform=False), name="$\\alpha$"))

# Have to add class and morphology columns at this stage otherwise it comes before these two in the LaTeX table
t.add_column(Column(ancs, name="Ancillary"))
t.add_column(Column(morphs, name="Morphology"))
t.add_column(Column(classes, name="Class"))

# Sort into order of Galactic l
t.sort(["GLON","GLAT"])
# Write out a copy for debugging, and for Vizier
t.write("concatenated_table.vot", format="votable", overwrite=True)

# Remove dead columns for LaTeX version
t.remove_column("GLON")
t.remove_column("GLAT")
t.remove_column("S200")
t.remove_column("err_S200")
t.remove_column("alpha")
t.remove_column("err_alpha")
     
t.write("table_for_my_snrs_paper.tex", format="latex", overwrite=True)
