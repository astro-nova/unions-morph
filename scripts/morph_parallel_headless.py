import warnings
warnings.filterwarnings('ignore')

import os
os.environ['PYTHONWARNINGS'] = 'ignore'

import numpy as np
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from functools import partial
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats

import sys
sys.path.append('/arc/home/esazonova/unions-morph')
from lib.io import download_files, make_cutout, parse_morph
from statmorph_lsst import SourceMorphology

import logging
import psutil
import time
from datetime import datetime

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('/arc/home/esazonova/unions-morph/logs/processing.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

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
        
    # print(f'{np.sum(good)} galaxies in the tile')
    logger.info(f'{np.sum(good)} galaxies in the tile')
    return cat, cat[good]


def process_galaxy(args):
    """Process a single galaxy from a tile"""
    idx, row, tilename, tile_f, weightmap_f, segmap_f = args
    try:
        
        # Make a cutout
        img, err, segmap, mask, psf, bgsd = make_cutout(row, tile_f, weightmap_f, segmap_f, r_frac=4)
        
        # Run statmorph
        isophotes = np.arange(22, 26.5, 0.5)
        pxscale = 0.1857  # arcsec/pixel
        fluxes = pxscale**2 * np.power(10, -(isophotes-30)/2.5)
        
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

        # Write result to file immediately
        out_file = '/arc/home/esazonova/unions-morph/catalogs/morph_new2.csv'
        with open(out_file, 'a') as f:
            # If filesize is 0 write header
            if f.tell() == 0:
                f.write(','.join(res.keys()) + '\n')
            f.write(','.join([str(v) for v in res.values()]) + '\n')

        logger.info(f"Done tile {tilename} galaxy {idx}")
    except Exception as e:
        logger.error(f"Error processing galaxy {idx} in tile {tilename}: {str(e)}")


def process_tile(tile):

    try:
        print(f'Processing tile {tile}')
        tilename = tile[-9:-2]

        # Download the data from arc to scratch
        download_files(tilename, tile=True, catalog=True, weightmap=True, segmap=True, photoz=False, star_galaxy=True)

        # Select galaxies
        cat, sample = select_sample(tilename, plot=False)

        # Open files and load data into memory
        tile_f = fits.open(f'/scratch/tile_{tilename}.fits')
        weightmap_f = fits.open(f'/scratch/wht_{tilename}.fits')
        segmap_f = fits.open(f'/scratch/seg_{tilename}.fits')
        
        # Prepare arguments for parallel galaxy processing
        galaxy_args = [
            (idx, row, tilename, tile_f, weightmap_f, segmap_f)
            for idx, row in sample.iterrows()
        ]
        
        # Process galaxies in parallel
        n_cores = max(1, cpu_count() - 1)
        with Pool(n_cores) as pool:
            results = pool.map(process_galaxy, galaxy_args)
        
        # Close files immediately
        tile_f.close()
        weightmap_f.close()
        segmap_f.close()

        # Record that this tile is done
        with open('/arc/home/esazonova/unions-morph/catalogs/processed_tiles_new.csv', 'a') as f:
            f.write(tilename + '\n')

        # Remove data from scratch
        os.remove(f'/scratch/tile_{tilename}.fits')
        os.remove(f'/scratch/wht_{tilename}.fits')
        os.remove(f'/scratch/seg_{tilename}.fits')
        os.remove(f'/scratch/cat_{tilename}.cat')
        os.remove(f'/scratch/sg_{tilename}.cat')
        logger.info(f"Completed processing tile {tile}")
    except:
        logger.info(f'Error processing tile {tile}')

    


if __name__ == '__main__':

    # Load tile list
    tile_df = pd.read_csv('/arc/home/esazonova/unions-morph/catalogs/tiles_r.csv')
    # done = pd.read_csv('/arc/home/esazonova/unions-morph/catalogs/processed_tiles_new.csv', names=['coords'])
    # done['tile'] = 'CFIS_LSB.' + np.char.mod('%07.3f', done.coords.values).astype(str) + '.r'
    # tile_df = tile_df[~tile_df.tile.isin(done.tile)]

    # Process each tile sequentially
    for i, tile in enumerate(tile_df.tile.values):
        process_tile(tile)
        if i > 0:
            break


