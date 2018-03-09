#!/usr/bin/env python

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

def reproject_snrs(snrs):
    padding = 1.4
    colors = ["white", "red", "green", "blue"]
    fitsdir = "/home/tash/data/MWA/GLEAM/GP/Week4/rgb/"

    for color in colors:
        print "Reprojecting "+color+" image"
        if not os.path.exists(color):
            os.makedirs(color)
        if color != "white":
            if not os.path.exists(color+"/rpj"):
                os.makedirs(color+"/rpj")

        fitsfile = fitsdir+color+"_MOL.fits"

        hdu = fits.open(fitsfile)
        w = wcs.WCS(hdu[0].header) 
        try:
            pix2deg = hdu[0].header["CDELT2"]
        except KeyError:
            pix2deg = hdu[0].header["CD2_2"]

        # Using the astropy cut-out method

        for snr in snrs:
            print "Reprojecting "+snr.name
    # No point doing the ones not in my search region
            l = snr.loc.galactic.l.value
            if (((l>180) and (l<240)) or (l>340) or (l<60)):

    # A bit of a sad hack, I'd rather just have one PSF file, but it's too much hassle
    # Week4 for GC SNR; Week2 for anticentre SNR
                if ((l>180) and (l<240)):
                    psf = fits.open(fitsdir+"Week2_"+color+"_lownoise_comp_psf.fits")
                else:
                    psf = fits.open(fitsdir+"Week4_"+color+"_lownoise_comp_psf.fits")
                name = snr.name+".fits"
# Set a minimum cutout size; otherwise I can't properly measure the remnant against the background
                if snr.min*60 > 20: 
                    framesize = u.Quantity(2*padding*snr.maj, u.deg)
                else:
                    framesize = u.Quantity(1, u.deg)
                cutout = Cutout2D(hdu[0].data, snr.loc, framesize, wcs=w)
            # Read these from the correct PSF image and then put them in the cutout
                wpsf = wcs.WCS(psf[0].header) 
                cp = np.squeeze(wpsf.wcs_world2pix([[snr.loc.ra.value,snr.loc.dec.value,1]],0))
                xp, yp = cp[0], cp[1]
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
                header_new["FREQ"] = hdu[0].header["FREQ"]
                new = fits.PrimaryHDU(cutout.data,header=header_new) #create new hdu
                newlist = fits.HDUList([new]) #create new hdulist
                newlist.writeto(color+"/"+name, overwrite = True)
                
            # Reproject the other colours
                if color != "white":
                    montage.mGetHdr("white/"+name,"temp.txt")
                    oldfile = color+"/"+name
                    newfile = color+"/rpj/"+name
                    montage.reproject(oldfile,newfile,header="temp.txt",exact_size=True)
# This copies ALL the fits keys across, including the beam values! So we need to replace those with the original values
                    oldhdr = fits.open(oldfile)[0].header
                    newhdu = fits.open(newfile)
                    newhdr = newhdu[0].header
                    for fitskey in ["BMAJ", "BMIN", "BPA", "FREQ"]:
                        newhdr[fitskey] = oldhdr[fitskey]
                        newhdu.writeto(newfile, overwrite = True)

if __name__ == "__main__":
    # Load the snrs fresh from the files
    snrs = read_snrs()
    reproject_snrs(snrs)
