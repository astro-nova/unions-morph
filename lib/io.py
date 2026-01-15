from vos import Client
import numpy as np
from skimage import transform as T
from astropy.stats import sigma_clipped_stats
from statmorph_lsst import _quantity_names
import os
import subprocess

def download_files(coords, weightmap=False, segmap=False, tile=False, catalog=False, 
    star_galaxy=False, photoz=False, path='/scratch/'):
    """ Download UNIONS files from VOSpace based on the provided coordinates. Optional args
    specify which files to download. `coords` should be in the XXX.YYY format following the 
    UNIONS tile naming convention. Downloads to scratch unless path is specified."""

    vosclient = Client()

    tilename = f'CFIS_LSB.{coords}.r'
    if tile:
        # download if doesn't exist already
        if not os.path.exists(f'{path}/tile_{coords}.fits'):
            vosclient.copy(f'vos:cfis/tiles_LSB_DR5/{tilename}.fits', 
                        f'{path}/tile_{coords}.fits')

    # Weightmap
    if weightmap:
        if not os.path.exists(f'{path}/wht_{coords}.fits'):
            vosclient.copy(f'vos:cfis/tiles_LSB_DR5/{tilename}.weight.fits.fz', 
                    f'{path}/wht_{coords}.fits.fz')
        
            # Decompress it
            subprocess.run(['funpack', f'{path}/wht_{coords}.fits.fz'], check=True)
            os.remove(f'{path}/wht_{coords}.fits.fz')

    # Catalog
    if catalog:
        if not os.path.exists(f'{path}/cat_{coords}.cat'):
            vosclient.copy(f'vos:cfis/tiles_DR5/CFIS.{coords}.r.cog.cat', 
                    f'{path}/cat_{coords}.cat')

    # Star-galaxy separation
    if star_galaxy:
        if not os.path.exists(f'{path}/sg_{coords}.cat'):
            vosclient.copy(f'vos:cfis/Processed_catalogues/StellarClass/stargal.cfis.r.dr5/tile.cats/CFIS.{coords}.r.sg.fits', 
                    f'{path}/sg_{coords}.cat')
    
    # Segmentation map
    if segmap:
        if not os.path.exists(f'{path}/seg_{coords}.fits'):
            vosclient.copy(f'vos:cfis/tiles_DR5/CFIS.{coords}.r.seg.fits.fz', 
                f'{path}/seg_{coords}.fits.fz')
            subprocess.run(['funpack', f'{path}/seg_{coords}.fits.fz'], check=True)
            os.remove(f'{path}/seg_{coords}.fits.fz')
        
    # Photoz catalog
    if photoz:
        if not os.path.exists(f'{path}/photoz_{coords}_ugriz.cat'):
            try:
                vosclient.copy(f'vos:cfis/gaap/UNIONS.{coords}_ugriz_photoz_ext.cat', 
                        f'{path}/photoz_{coords}_ugriz.cat')
            except HTTPError:
                vosclient.copy(f'vos:cfis/gaap/UNIONS.{coords}_ugri_photoz_ext.cat', 
                        f'{path}/photoz_{coords}_ugriz.cat')
            except:
                print('No photo-z catalog available yet')

def make_cutout(galaxy, tile, weightmap, segmap, cutout_min=20, r_frac=2):
    """ Make cutouts of a galaxy from the tile, weightmap and segmentation map.
    `galaxy` is a row from the catalog dataframe. `tile`, `weightmap` and `segmap` are  
    the opened fits files. `cutout_min` is the minimum size of the cutout in pixels.
    `r_frac` is the factor multiplied by FLUX_RADIUS to define the cutout size. """

    # Cutout size
    xc, yc = int(galaxy.X_IMAGE+0.5), int(galaxy.Y_IMAGE+0.5)
    axis_ratio = np.max([0.2, galaxy.Q])
    size = int(galaxy.FLUX_RADIUS*r_frac / axis_ratio )
    size = np.max([size, cutout_min])

    # If the cutout goes beyond the image edges, adjust the slices.
    # Compact formulation. 
    y_start, y_end = max(0, yc-size), min(tile[0].data.shape[0], yc+size)
    x_start, x_end = max(0, xc-size), min(tile[0].data.shape[1], xc+size)
    slices = slice(y_start, y_end), slice(x_start, x_end)

    # Make cutouts 
    img = tile[0].data[slices]
    segmap = segmap[1].data[slices]
    err = weightmap[1].data[slices]
    mask = err==0

    # Only select the source in the segmap
    xc = img.shape[1]//2
    yc = img.shape[0]//2
    source = segmap[yc, xc]
    segmap[(segmap > 0) & (segmap != source)] = 2
    segmap[segmap == source] = 1

    # Estimate and subtract a flat background
    bgmean, bgmed, bgsd = sigma_clipped_stats(img, mask=mask|(segmap>0))
    img -= bgmed

    # Load the empirical PSF
    psf_1arcsec = np.load(f'/arc/home/esazonova/unions-morph/data/psf_1arcsec.npy')
    fwhm = galaxy.PREDIQ
    # Generate the PSF with that FWHM from the psfex fit
    psf = T.rescale(psf_1arcsec, (fwhm/1))
    psf = psf/np.sum(psf)
    
    return img, err, segmap, mask, psf, bgsd

def parse_morph(out_dict, morph):

    qs = _quantity_names
    qs += ['flag','flag_sersic']
    isophotes = np.arange(22, 26.5, 0.5)
    for q in qs:
        if q == 'isophote_asymmetry':
            aisos = morph.__getattribute__(q)
            for sblim, aiso in zip(isophotes, aisos):
                out_dict[f'aiso_{sblim:0.1f}'] = aiso
        else:
            out_dict[q] = morph.__getattribute__(q)
    return out_dict