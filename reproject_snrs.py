#!/usr/bin/env python

__author__ = "Natasha Hurley-Walker"
__date__ = "10/03/2018"

import os
import sys
from astropy.io import fits
from astropy import wcs
from astropy.nddata import Cutout2D
from astropy import units as u

import numpy as np

# Regridding
import montage_wrapper as montage

from defs import Coords
from defs import SNR
from defs import read_snrs

# For each remnant, cut out the region around it in the white fits file
# Then montage the red, green, and blue fits files to match the white

def reproject_snrs(snrs, clobber):
    padding = 4.0
#    colors = ["072-080MHz", "080-088MHz", "088-095MHz", "095-103MHz", "103-111MHz", "111-118MHz", "118-126MHz", "126-134MHz", "139-147MHz", "147-154MHz", "154-162MHz", "162-170MHz", "170-177MHz", "177-185MHz", "185-193MHz", "193-200MHz", "200-208MHz", "208-216MHz", "216-223MHz", "223-231MHz"]
    colors = ["white", "red", "green", "blue", "072-080MHz", "080-088MHz", "088-095MHz", "095-103MHz", "103-111MHz", "111-118MHz", "118-126MHz", "126-134MHz", "139-147MHz", "147-154MHz", "154-162MHz", "162-170MHz", "170-177MHz", "177-185MHz", "185-193MHz", "193-200MHz", "200-208MHz", "208-216MHz", "216-223MHz", "223-231MHz"]
    fitsdir = "/home/tash/data/MWA/GLEAM/GP/"

    for color in colors:
        print "Reprojecting "+color+" image"
        if not os.path.exists(color):
            os.makedirs(color)
        if color != "white":
            if not os.path.exists(color+"/rpj"):
                os.makedirs(color+"/rpj")

#        fitsfile = fitsdir+color+"_MOL.fits"
#
#        hdu = fits.open(fitsfile)
#        w = wcs.WCS(hdu[0].header) 
#        try:
#            pix2deg = hdu[0].header["CDELT2"]
#        except KeyError:
#            pix2deg = hdu[0].header["CD2_2"]

        # Using the astropy cut-out method

        for snr in snrs:
            print "Reprojecting "+snr.name
    # No point doing the ones not in my search region
            l = snr.loc.galactic.l.value
            if (((l>180) and (l<240)) or (l>300) or (l<60)):
                name = snr.name+".fits"
                if clobber or not os.path.exists(color+"/"+name) or not os.path.exists(color+"/rpj/"+name):
        # Week4 for GC SNR; Week2 for anticentre SNR
                    if ((l>180) and (l<240)):
                        psf = fits.open(fitsdir+"Week2/Week2_"+color+"_lownoise_comp_psf.fits")
                        hdu = fits.open(fitsdir+"Week2/Week2_"+color+"_lownoise_ddmod_rescaled.fits")
# HACK
#                        orig_hdu = fits.open(fitsdir+"Week2_"+color+"_lownoise_ddmod_rescaled.fits")
#                        hdu = fits.open("/home/tash/data/MWA/GLEAM/allmosaics/SNR_G189.6+3.3/"+color+"/"+snr.name+".fits")
                    else:
                        psf = fits.open(fitsdir+"Week4/Week4_"+color+"_lownoise_comp_psf.fits")
                        hdu = fits.open(fitsdir+"Week4/Week4_"+color+"_lownoise_ddmod_rescaled.fits")
                    try:
                        pix2deg = hdu[0].header["CDELT2"]
                    except KeyError:
                        pix2deg = hdu[0].header["CD2_2"]
                    w = wcs.WCS(hdu[0].header) 
    # Set a minimum cutout size; otherwise I can't properly measure the remnant against the background
                    if snr.min*60 > 20: 
                        framesize = u.Quantity(2*padding*snr.maj, u.deg)
                    else:
                        framesize = u.Quantity(1, u.deg)
   # Set a maximum cutout size to avoid getting weird results for some snr
                    if framesize.value > 4.0:
                        framesize = u.Quantity(4, u.deg)
                       
                    print framesize
                    cutout = Cutout2D(hdu[0].data, snr.loc.fk5, framesize, wcs=w)
                # Read these from the correct PSF image and then put them in the cutout
                    wpsf = wcs.WCS(psf[0].header) 
                    cp = np.squeeze(wpsf.wcs_world2pix([[snr.loc.fk5.ra.value,snr.loc.fk5.dec.value,1]],0))
                    xp, yp = int(cp[0]), int(cp[1])
                    bmaj = psf[0].data[0,yp,xp]
                    bmin = psf[0].data[1,yp,xp]
                    bpa = psf[0].data[2,yp,xp]
            #    beamvolume = (1.1331 * bmaj * bmin) # gaussian beam conversion
                    header_new = cutout.wcs.to_header()
                # Edit the header so that the CD values are copied from the PC values -- makes it DS9-readable
                    header_new["CD1_1"] = header_new["PC1_1"]
                    header_new["CD2_2"] = header_new["PC2_2"]
                    header_new["BMAJ"] = bmaj
                    header_new["BMIN"] = bmin
                    header_new["BPA"] = bpa
#HACK
                    header_new["FREQ"] = hdu[0].header["FREQ"]
#                    header_new["FREQ"] = orig_hdu[0].header["FREQ"]
                    new = fits.PrimaryHDU(cutout.data,header=header_new) #create new hdu
                    newlist = fits.HDUList([new]) #create new hdulist
                    try:
                        newlist.writeto(color+"/"+name, overwrite = True)
# For some reason I get:
# NameError: global name 'VerifyError' is not defined
                    except:
                        print "Invalid fits keys for {0} at {1}".format(snr.name,color)
                    
                    print snr.name, "here"
                # Reproject the other colours
                    if color != "white":
                        montage.mGetHdr("white/"+name,"temp.txt")
                        oldfile = color+"/"+name
                        newfile = color+"/rpj/"+name
                        try:
                            montage.reproject(oldfile,newfile,header="temp.txt",exact_size=True)
                            oldhdr = fits.open(oldfile)[0].header
                            newhdu = fits.open(newfile)
                            newhdr = newhdu[0].header
# This copies ALL the fits keys across, including the beam values! So we need to replace those with the original values
                            for fitskey in ["BMAJ", "BMIN", "BPA", "FREQ"]:
                                newhdr[fitskey] = oldhdr[fitskey]
                            newhdu.writeto(newfile, overwrite = True)
# For some reason I get:
# NameError: global name 'MontageError' is not defined
#                        except MontageError:
                        except:
                            print "Montage reprojection failed for {0} at {1}".format(snr.name,color)

if __name__ == "__main__":
    # Load the snrs fresh from the files
    snrs = read_snrs()
    clobber = True
    reproject_snrs(snrs, clobber)
