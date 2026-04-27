import warnings
warnings.filterwarnings('ignore')

import os
os.environ['PYTHONWARNINGS'] = 'ignore'
from multiprocessing import Pool, cpu_count
from concurrent.futures import ThreadPoolExecutor
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
import fcntl
import traceback

from functools import partial
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats

import sys

# Resolve paths from the environment so the same script runs on CANFAR, Nibi,
# and laptops. PROJECT_DIR falls back to the repo root (scripts/'s parent).
PROJECT_DIR = os.environ.get(
    'PROJECT_DIR',
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)
WORK_DIR = os.environ.get('WORK_DIR', '/scratch')
OUT_DIR = os.environ.get('OUT_DIR', os.path.join(PROJECT_DIR, 'catalogs'))
LOG_DIR = os.environ.get('LOG_DIR', os.path.join(PROJECT_DIR, 'logs'))
os.makedirs(WORK_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

sys.path.insert(0, PROJECT_DIR)
from lib.io import make_cutout, parse_morph, download_files
from statmorph_lsst import SourceMorphology

import logging
import psutil
import signal
import time
from datetime import datetime


GALAXY_TIMEOUT_SECONDS = 120


class _GalaxyTimeout(Exception):
    pass


def _galaxy_alarm(signum, frame):
    raise _GalaxyTimeout(f"statmorph exceeded {GALAXY_TIMEOUT_SECONDS}s budget")

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(os.path.join(LOG_DIR, 'processing.log')),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)



def select_sample(tile, plot=False):

    # Load the source extractor and star-galaxy catalogs
    cat = Table.read(f'{WORK_DIR}/cat_{tile}.cat', format='ascii').to_pandas()
    sg = Table.read(f'{WORK_DIR}/sg_{tile}.cat').to_pandas()[['M1','M2','M3','s21','s31']]
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
    idx, row, tilename, psf = args
    
    # Use variables to track what needs cleanup
    img = err = segmap = mask = morph = None
    
    try:
        # Use context managers to ensure files are always closed
        with fits.open(f'{WORK_DIR}/tile_{tilename}.fits', memmap=True, lazy_load_hdus=True) as tile_f, \
             fits.open(f'{WORK_DIR}/wht_{tilename}.fits', memmap=True, lazy_load_hdus=True) as weightmap_f, \
             fits.open(f'{WORK_DIR}/seg_{tilename}.fits', memmap=True, lazy_load_hdus=True) as segmap_f:

            # Make a cutout
            img, err, segmap, mask, psf, bgsd = make_cutout(row, tile_f, weightmap_f, segmap_f, r_frac=4, psf=psf)
        
        # Run statmorph (files are now closed)
        isophotes = np.arange(22, 26.5, 0.5)
        pxscale = 0.1857  # arcsec/pixel
        fluxes = pxscale**2 * np.power(10, -(isophotes-30)/2.5)
        
        signal.signal(signal.SIGALRM, _galaxy_alarm)
        signal.alarm(GALAXY_TIMEOUT_SECONDS)
        try:
            morph = SourceMorphology(
                img, segmap, label=1, weightmap=err, mask=mask, psf=psf,
                interpolate_mask=False, asymmetry_isophotes=fluxes,
                sersic_model_args={'bounds' : {'n' : (0.1, 6)}}
            )
        finally:
            signal.alarm(0)
        
        # Parse output
        res = {
            'tile' : tilename, 'idx' : row.name, 'ra': row['ALPHA_J2000'], 
            'dec' : row['DELTA_J2000'], 'fwhm' : row.PREDIQ}
        res = parse_morph(res, morph)

        # Write result to file immediately with locking to prevent race conditions
        out_file = f'{WORK_DIR}/morph_batch_{tilename}.csv'
        with open(out_file, 'a') as f:
            # Acquire exclusive lock before writing
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            try:
                # Check if file is empty and write header if needed
                f.seek(0, 2)  # Seek to end
                if f.tell() == 0:
                    f.write(','.join(res.keys()) + '\n')
                f.write(','.join([str(v) for v in res.values()]) + '\n')
                f.flush()  # Ensure data is written
            finally:
                # Release lock
                fcntl.flock(f.fileno(), fcntl.LOCK_UN)

        # Explicitly clean up large objects to free memory
        del img, err, segmap, mask, psf, morph
        
        return True  # Success indicator

    except Exception as e:
        # Clean up any large objects that were created before the error
        # This prevents BrokenPipeError from trying to pickle large objects in the exception
        try:
            if img is not None: del img
            if err is not None: del err
            if segmap is not None: del segmap
            if mask is not None: del mask
            if psf is not None: del psf
            if morph is not None: del morph
        except:
            pass
        
        # Log error with minimal info to avoid pickling issues
        logger.error(f"Error processing galaxy {idx} in tile {tilename}: {str(e)[:200]}s")

        
        # Return False to indicate failure (don't let exception propagate to multiprocessing)
        return False


def preprocess_tile(tile, psf):
    """Download a tile and build the per-galaxy work list. Runs on the
    main/downloader core; overlaps with the previous tile's galaxy workers."""
    tilename = tile[-9:-2]
    logger.info(f"Preprocessing tile {tile}")

    download_files(tilename, tile=True, catalog=True, weightmap=True,
                   segmap=True, photoz=False, star_galaxy=True, path=WORK_DIR)

    cat, sample = select_sample(tilename, plot=False)
    galaxy_args = [(idx, row, tilename, psf) for idx, row in sample.iterrows()]
    del cat, sample
    return {'tile': tile, 'tilename': tilename, 'galaxy_args': galaxy_args}


def finalize_tile(info):
    """Record completion, clean scratch files, and write the final fits.
    Runs on the main core after the galaxy workers for this tile finish."""
    tile = info['tile']
    tilename = info['tilename']

    with open(os.path.join(OUT_DIR, 'processed_tiles_new2.csv'), 'a') as f:
        f.write(tile + '\n')

    for p in (f'{WORK_DIR}/tile_{tilename}.fits',
              f'{WORK_DIR}/wht_{tilename}.fits',
              f'{WORK_DIR}/seg_{tilename}.fits',
              f'{WORK_DIR}/cat_{tilename}.cat',
              f'{WORK_DIR}/sg_{tilename}.cat'):
        if os.path.exists(p):
            os.remove(p)

    morph_csv = f'{WORK_DIR}/morph_batch_{tilename}.csv'
    if os.path.exists(morph_csv):
        morph_cat = pd.read_csv(morph_csv)
        morph_table = Table.from_pandas(morph_cat)
        morph_table.write(os.path.join(OUT_DIR, f'morphology/morph_tile_{tilename}.fits'),
                          overwrite=True)

    logger.info(f"Completed processing tile {tile}")


def mark_error(tilename):
    with open(os.path.join(OUT_DIR, 'error_tiles.csv'), 'a') as f:
        f.write(tilename + '\n')




if __name__ == '__main__':


    # Parse input args: imin, imax to index tile df using argparse
    parser = argparse.ArgumentParser(description='Process tiles for morphology analysis.')
    parser.add_argument('--imin', type=int, required=True,  help='Minimum index of tiles to process')
    parser.add_argument('--imax', type=int, required=True,  help='Maximum index of tiles to process')
    args = parser.parse_args()

    # Load tile list, slice by index, then drop tiles already in the done log.
    tile_df = pd.read_csv(os.path.join(PROJECT_DIR, 'catalogs', 'tiles_todo.csv'))
    tile_df = tile_df.iloc[args.imin:args.imax]

    done_path = os.path.join(OUT_DIR, 'processed_tiles_new2.csv')
    if os.path.exists(done_path):
        done = pd.read_csv(done_path, names=['tile'])
        n_before = len(tile_df)
        tile_df = tile_df[~tile_df.tile.isin(done.tile)]
        logger.info(
            f"Skipping {n_before - len(tile_df)} already-processed tiles "
            f"(from {done_path})"
        )

    tiles = list(tile_df.tile.values)

    # Load the PSF once and share it across tiles (it's static).
    psf = np.load(os.path.join(PROJECT_DIR, 'data', 'psf_1arcsec.npy'))

    # Respect the SLURM cpus-per-task allocation (cpu_count() reports the whole
    # node, not our share). Reserve one core for the downloader/main process;
    # the rest run galaxy workers.
    n_cores = int(os.environ.get(
        'SLURM_CPUS_PER_TASK',
        len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else cpu_count()
    ))
    n_galaxy_cores = max(1, n_cores - 1)

    # Downloads are network-bound (~9 min/tile, mostly waiting), so multiple
    # can run concurrently without touching CPU. Galaxy processing on 89 cores
    # finishes a tile in ~2-3 min, so we need ~3-4 downloads in flight to keep
    # workers saturated. Configurable via N_DOWNLOADERS env var.
    n_downloaders = int(os.environ.get('N_DOWNLOADERS', 4))
    logger.info(
        f"Pipeline: {n_galaxy_cores} galaxy workers + {n_downloaders} "
        f"concurrent downloaders, {len(tiles)} tiles to process"
    )

    # Fork the persistent worker pool now, before any tile data is loaded, so
    # workers don't inherit preprocessing state. The same pool is reused across
    # tiles: we submit each tile's galaxies as one async batch.
    pool = Pool(n_galaxy_cores)
    try:
        # Kick off all downloads up front. The ThreadPoolExecutor throttles
        # to n_downloaders concurrent; the rest queue and start as slots free.
        with ThreadPoolExecutor(max_workers=n_downloaders) as dl_ex:
            dl_futures = [dl_ex.submit(preprocess_tile, t, psf) for t in tiles]

            # Consume in submission order. fut.result() blocks on THIS tile's
            # download, but other downloads keep running in the background.
            for i, fut in enumerate(dl_futures):
                try:
                    info = fut.result()
                except Exception:
                    logger.error(
                        f"Preprocess failed for tile {tiles[i]}:\n{traceback.format_exc()}"
                    )
                    mark_error(tiles[i][-9:-2])
                    continue

                # Dispatch this tile's galaxies to the 89-worker pool.
                try:
                    if info['galaxy_args']:
                        map_result = pool.map_async(
                            process_galaxy, info['galaxy_args'], chunksize=1
                        )
                        map_result.wait()
                    finalize_tile(info)
                except Exception:
                    logger.error(
                        f"Processing failed for tile {info['tile']}:\n{traceback.format_exc()}"
                    )
                    mark_error(info['tilename'])
    finally:
        pool.close()
        pool.join()