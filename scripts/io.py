from vos import Client
import sys
sys.path.append('/home/u1/r/rmjarvis/py/lib/python3.8/site-packages/')

def download_files(coords, weightmap=False, segmap=False, tile=False, catalog=False, 
    star_galaxy=False, photoz=False, path='/scratch/'):
    """ Download UNIONS files from VOSpace based on the provided coordinates. Optional args
    specify which files to download. `coords` should be in the XXX.YYY format following the 
    UNIONS tile naming convention. Downloads to scratch unless path is specified."""

    vosclient = Client()

    tilename = f'CFIS_LSB.{tile}.r'
    if tile:
        vosclient.copy(f'vos:cfis/tiles_LSB_DR5/{tilename}.fits', 
                       f'{path}/tile_{coords}.fits')

    # Weightmap
    if weightmap:
        vosclient.copy(f'vos:cfis/tiles_LSB_DR5/{tilename}.weight.fits.fz', 
                   f'{path}/wht_{coords}.fits.fz')
        # Decompress it
        os.system(f'funpack {path}/wht_{coords}.fits.fz')

    # Catalog
    if catalog:
        vosclient.copy(f'vos:cfis/tiles_DR5/CFIS.{coords}.r.cog.cat', 
                   f'{path}/cat_{coords}.cat')

    # Star-galaxy separation
    if star_galaxy:
        vosclient.copy(f'vos:cfis/Processed_catalogues/StellarClass/stargal.cfis.r.dr5/tile.cats/CFIS.{tile}.r.sg.fits', 
                   f'{path}/sg_{coords}.cat')
    
    # Segmentation map
    if segmap:
        vosclient.copy(f'vos:cfis/tiles_DR5/CFIS.{coords}.r.seg.fits.fz', 
                   f'{path}/seg_{coords}.fits.fz')
        os.system(f'funpack {path}/seg_{coords}.fits.fz')
    
    # Photoz catalog
    if photoz:
        try:
            vosclient.copy(f'vos:cfis/gaap/UNIONS.{coords}_ugriz_photoz_ext.cat', 
                       f'{path}/photoz_{coords}_ugriz.cat')
        except HTTPError:
            vosclient.copy(f'vos:cfis/gaap/UNIONS.{coords}_ugri_photoz_ext.cat', 
                       f'{path}/photoz_{coords}_ugriz.cat')
        except:
            print('No photo-z catalog available yet')