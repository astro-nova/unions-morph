import os
import numpy as np
import pandas as pd
from tqdm.notebook import tqdm
from joblib import parallel, delayed
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats

import sys
sys.path.append('..')
from lib.io import download_files, make_cutout, parse_morph
from statmorph_lsst import SourceMorphology

# Suppress warnings from photutils, astropy or statmorph_lsst
import warnings
warnings.filterwarnings('ignore')

def select_sample(tile, plot=False):

    # Load the source extractor and star-galaxy catalogs
    cat = Table.read(f'/scratch/cat_{tile}.cat', format='ascii').to_pandas()
    sg = Table.read(f'/scratch/sg_{tile}.cat').to_pandas()[['M1','M2','M3','s21','s31']]
    cat = pd.merge(cat, sg, left_index=True, right_index=True)

    # Calculate surface brightness and axis ratio proxies
    cat['Q'] = cat['B_WORLD']/cat['A_WORLD']
    
    # Apply criteria above
    good = ((cat.FLUX_RADIUS >= 4) &
            (cat.FLAGS < 17) &
            (cat.Q >= 0.05) &
            (cat.MAG_COG <= 27) &
            (cat.MAG_COG >= 14) &
            (cat.M1-cat.M2 < 1.5) &
            (cat.M1-cat.M2 > 0.5) 
           )
    
    # Also from the star-galaxy catalog
    stars = (cat.MAG_COG < 23.5) & (cat.M1-cat.M2 >= 0.63) & (cat.M1-cat.M2 < 0.75)
    stars = (cat.s21 < 3) | (cat.s31 < 3)
    
    # Plot the selected sample
    if plot:
        fig, ax = plt.subplots(1,1, figsize=(8,5))
        plt.scatter(cat.MAG_COG, cat.M1-cat.M2, s=1, alpha=0.1, c='k')
        plt.scatter(cat[good].MAG_COG, cat[good].M1-cat[good].M2, s=1, alpha=0.5, c='b')
        plt.scatter(cat[stars].MAG_COG, cat[stars].M1-cat[stars].M2, s=1, alpha=0.5, c='r')
        plt.xlim(14, 28)
        plt.ylim(-0.5, 1.5)
        plt.colorbar()

    good = good & ~stars
        
    print(f'{np.sum(good)} galaxies in the tile')
    return cat, cat[good]

def process_tile(tilerow):

    tilename = tilerow.tile[-9:-2]

    # Download the data from arc to scratch
    download_files(tilename, tile=True, catalog=True, weightmap=True, segmap=True, photoz=False, star_galaxy=True)

    # Select galaxies
    cat, sample = select_sample(tilename, plot=False)

    # Open files
    tile_f = fits.open(f'/scratch/tile_{tilename}.fits')
    weightmap_f = fits.open(f'/scratch/wht_{tilename}.fits')
    segmap_f = fits.open(f'/scratch/seg_{tilename}.fits')

    # For each galaxy, make a cutout and run statmorph
    isophotes = np.arange(22, 26.5, 0.5)
    pxscale = 0.1857  # arcsec/pixel
    fluxes = pxscale**2 * np.power(10, -(isophotes-30)/2.5)

    for idx, row in sample.iterrows():
        try:
            # Make a cutout
            img, err, segmap, mask, psf, bgsd = make_cutout(row, tile_f, weightmap_f, segmap_f, r_frac=4)
            # Run statmorph
            morph = SourceMorphology(
                img, segmap, label=1, weightmap=err, mask=mask, psf=psf, 
                interpolate_mask=False, asymmetry_isophotes=fluxes,
                sersic_model_args={'bounds' : {'n' : (0.1, 6)}}
            )
            # Parse output
            res = {
                'tile' : tilename, 'idx' : row.name, 'ra': row['ALPHA_J2000'], 
                'dec' : row['DELTA_J2000'], 'fwhm' : row.PREDIQ}
            res = parse_morph(res, morph)

            # Write the output
            out_file = f'../catalogs/morph.csv'
            with open(out_file, 'a') as f:
                # If filesize is 0 write header
                if f.tell() == 0:
                    f.write(','.join(res.keys()) + '\n')
                f.write(','.join([str(v) for v in res.values()]) + '\n')
        except:
            continue

    # Close files
    tile_f.close()
    weightmap_f.close()
    segmap_f.close()

    # Record that this tile is done
    with open('../catalogs/processed_tiles.csv', 'a') as f:
        f.write(tilename + '\n')

    # Remove data from scratch
    os.remove(f'/scratch/tile_{tilename}.fits')
    os.remove(f'/scratch/wht_{tilename}.fits')
    os.remove(f'/scratch/seg_{tilename}.fits')
    os.remove(f'/scratch/cat_{tilename}.cat')
    os.remove(f'/scratch/sg_{tilename}.cat')
    


if __name__ == '__main__':

    tile_df = pd.read_csv('../catalogs/tiles_r.csv')

    # Start delayed joblib run
    # parallel(n_jobs=16)(
    #     delayed(process_tile)(row) for idx, row in tqdm(tile_df.iterrows(), total=len(tile_df))
    # )
    process_tile(tile_df.iloc[0])


