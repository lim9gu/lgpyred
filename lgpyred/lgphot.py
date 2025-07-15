### Photometry

import time
import os, glob
import subprocess

import numpy as np
from numpy import median

import astropy.units as u
from astropy.io import ascii, fits
from astropy.wcs import WCS
from astropy.table import vstack, Table
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize)
from astropy.visualization import ZScaleInterval, LinearStretch

from importlib.resources import files
import importlib.resources as pkg_resources
import lgpyred.imsng.imsngphot as iph 
from lgpyred.reduction.hdrcheck import wcscenter
import lgpyred.data

from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)


def PlotPhotRes(inim, band, ref, mykey, myerrkey, totals, alives, exiles, ZP, zpstd, fwhm, FWHM, FWHM_std):
    #PlotPhotRes(inim, band, ref, mykey, myerrkey, totals, alives, exiles, ZP_dict[mykey], ZPer_dict[mykey], fwhm, FWHM, FWHM_std)

    width = 1.5
    rc('axes', linewidth=width)
    plt.rcParams["font.family"] = "monospace"
    fig  = plt.figure(0, figsize=(14, 12))
    grid = GridSpec(ncols=2, nrows=2)

    ax1  = fig.add_subplot(grid[0, 0])  
    ax1.tick_params(axis='both', which='minor', direction='in', labelsize=13, labelbottom=True, width=width, length=4)
    ax1.tick_params(axis='both', which='major', direction='in', labelsize=13, labelbottom=True, width=width, length=8)
    ax1.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax1.xaxis.set_major_locator(MultipleLocator(1))
    ax1.xaxis.set_ticks_position('both')
    ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax1.yaxis.set_major_locator(MultipleLocator(0.2))
    ax1.yaxis.set_ticks_position('both') 
    ax1.errorbar(totals[band], totals['zp'], yerr=totals['zper'], color='black', fmt='s', mfc='white', capsize=3) # total
    ax1.errorbar(alives[band], alives['zp'], yerr=alives['zper'], color='black', fmt='s', mfc='dodgerblue', capsize=3)
    ax1.errorbar(exiles[band], exiles['zp'], yerr=exiles['zper'], color='black', fmt='s', mfc='orangered', capsize=3)
    ax1.hlines(ZP, np.min(totals[band])-2, np.max(totals[band])+2, linestyle='solid', color='black', linewidth=1.5 )
    ax1.hlines(ZP+1*zpstd, np.min(totals[band])-2, np.max(totals[band])+2, linestyle='dashed', color='black', linewidth=1 )
    ax1.hlines(ZP-1*zpstd, np.min(totals[band])-2, np.max(totals[band])+2, linestyle='dashed', color='black', linewidth=1 )
    ax1.text( 0.9, 0.9, f'ZP={ZP:.3f}+/-{zpstd:.3f}', ha='right', va='top', transform=ax1.transAxes, fontsize=12)
    ax1.text( 0.9, 0.83, f'# of good={len(alives)}', ha='right', va='top', transform=ax1.transAxes, fontsize=12)
    ax1.text( 0.9, 0.76, f'{mykey}', ha='right', va='top', transform=ax1.transAxes, fontsize=12)

    ax1.set_xlabel(f'{ref} {band} (AB mag)', {'color':'black', 'fontsize':15})
    ax1.set_ylabel(f'ZP (mag)', {'color':'black', 'fontsize':15})
    ax1.set_xlim(np.min(totals[band])-1,  np.max(totals[band])+1)
    ax1.set_ylim(ZP-0.5, ZP+0.5)
    ax1.set_title("ZP measurement", fontsize=16)

    ax2  = fig.add_subplot(grid[0, 1]) 
    ax2.tick_params(axis='both', which='minor', direction='in', labelsize=13, labelbottom=True, width=width, length=4)
    ax2.tick_params(axis='both', which='major', direction='in', labelsize=13, labelbottom=True, width=width, length=8)  
    ax2.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax2.xaxis.set_major_locator(MultipleLocator(1))
    ax2.xaxis.set_ticks_position('both')
    ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax2.yaxis.set_major_locator(MultipleLocator(0.2))
    ax2.yaxis.set_ticks_position('both') 
    ax2.errorbar(totals[band], totals[band] - (ZP + totals['MAG_'+mykey]), yerr=np.sqrt(totals[band+'err']**2 + totals['zper']**2 + totals[myerrkey]**2), color='black', fmt='s', mfc='white', capsize=3)
    ax2.errorbar(alives[band], alives[band] - (ZP + alives['MAG_'+mykey]), yerr=np.sqrt(alives[band+'err']**2 + alives['zper']**2 + alives[myerrkey]**2), color='black', fmt='s', mfc='dodgerblue', capsize=3)
    ax2.errorbar(exiles[band], exiles[band] - (ZP + exiles['MAG_'+mykey]), yerr=np.sqrt(exiles[band+'err']**2 + exiles['zper']**2 + exiles[myerrkey]**2), color='black', fmt='s', mfc='orangered', capsize=3)
    ax2.hlines(0, np.min(totals[band])-2, np.max(totals[band])+2, linestyle='solid', color='black', linewidth=1.5 )
    ax2.hlines(0+1*zpstd, np.min(totals[band])-2, np.max(totals[band])+2, linestyle='dashed', color='black', linewidth=1 )
    ax2.hlines(0-1*zpstd, np.min(totals[band])-2, np.max(totals[band])+2, linestyle='dashed', color='black', linewidth=1 )
    ax2.set_xlabel(f'{ref} {band} (AB mag)', {'color':'black', 'fontsize':15})
    ax2.set_ylabel(f'{ref}-Mag. (AB mag)', {'color':'black', 'fontsize':15})
    ax2.set_xlim(np.min(totals[band])-1,  np.max(totals[band])+1)
    ax2.set_ylim(0-0.5, 0+0.5)
    ax2.set_title("ZP checkplot", fontsize=16)

    data, hdr    = fits.getdata(inim, header=True)
    wcs          = WCS(hdr)
    norm_zscale  = ImageNormalize(data, interval=ZScaleInterval(), stretch=LinearStretch())
    ax3  = fig.add_subplot(grid[1, 0], projection=wcs) 
    im   = ax3.imshow(data, cmap='gray', origin='lower', norm=norm_zscale)
    for x, y, n in zip(totals['X_IMAGE'], totals['Y_IMAGE'], totals['refnum']) :
        circ = Circle((x, y), 25, color='gold', fill=None, linewidth=0.5)
        ax3.add_patch(circ)
        ax3.text(x, y+50, str(n), color='gold', fontsize=10)
    for x, y, n in zip(alives['X_IMAGE'], alives['Y_IMAGE'], alives['refnum']) :
        circ = Circle((x, y), 25, color='blue', fill=None, linewidth=0.5)
        ax3.add_patch(circ)
        ax3.text(x, y+50, str(n), color='blue', fontsize=10)    
    for x, y, n in zip(exiles['X_IMAGE'], exiles['Y_IMAGE'], exiles['refnum']) :
        circ = Circle((x, y), 25, color='red', fill=None, linewidth=0.5)
        ax3.add_patch(circ)
        ax3.text(x, y+50, str(n), color='red', fontsize=10)  
    ax3.coords[0].set_axislabel("Right Ascension", fontsize=16)
    ax3.coords[1].set_axislabel("Declination", fontsize=16)
    ax3.set_title("Selected ref. stars", fontsize=16)

    ax4  = fig.add_subplot(grid[1, 1]) 
    ax4.tick_params(axis='both', which='minor', direction='in', labelsize=13, labelbottom=True, width=width, length=4)
    ax4.tick_params(axis='both', which='major', direction='in', labelsize=13, labelbottom=True, width=width, length=8)  
 
    h = ax4.hist(fwhm[fwhm < 7], bins=30, histtype='step', linewidth=1.5, color='black') 
    ax4.vlines(FWHM, 0, np.max(h[0]+1), linestyle='solid', color='dodgerblue', linewidth=1.5 )
    ax4.vlines(FWHM+1*FWHM_std, 0, np.max(h[0]+1), linestyle='dashed', color='dodgerblue', linewidth=1.5 )
    ax4.vlines(FWHM-1*FWHM_std, 0, np.max(h[0]+1), linestyle='dashed', color='dodgerblue', linewidth=1.5 )
    ax4.set_xlim(0,  7)
    ax4.set_ylim(0,  np.max(h[0]+1))
    ax4.text( 0.9, 0.9, f'FWHM={FWHM:.2f}+/-{FWHM_std:.2f}(")', ha='right', va='top', transform=ax4.transAxes, fontsize=12)
    ax4.set_title("FWHM histogram", fontsize=16)
    ax4.set_xlabel('FWHM (")', {'color':'black', 'fontsize':15})
    ax4.set_ylabel('Number', {'color':'black', 'fontsize':15})
    outname = inim[:-5]+'.png'

    fig.suptitle(outname, fontsize=18)

    fig.savefig(outname, dpi=200, facecolor='w', edgecolor='w', orientation='portrait', format=None, transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)
    #plt.savefig(outname, dpi=200, facecolor='w', edgecolor='w', orientation='portrait', format=None, transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)
    print(outname+' is created.')
    os.system('mkdir fig')
    #os.makedirs('fig', exist_ok=True)
    os.system(f'mv {outname} ./fig/')
    plt.close()

def SE(inim, outcat_name, pixscale, path=None, gain=1, FWHM=1.2, minarea=5, det_thresh=5, deb_nthresh=32, deb_mincont=0.01, backsize=128, backfiltersize=5, backphoto_type='LOCAL', backphoto_thick=24, aperture=None):

    if path is None:
        path = pkg_resources.resource_filename('lgpyred', 'photconf')

    basepath = files('lgpyred.photconf')
    sexconf  = basepath.joinpath('default.sex')
    sexnnw   = basepath.joinpath('default.nnw')
    sexconv  = basepath.joinpath('default.conv')
    sexparam = basepath.joinpath('default.param')

    checkimage_type = 'BACKGROUND,SEGMENTATION,APERTURES'
    checkimage_name = inim[:-5]+'.chbkg.fits,'+inim[:-5]+'.chseg.fits,'+inim[:-5]+'.chaper.fits'

    sexcom = f"sex {inim} -c {sexconf} -CATALOG_NAME {outcat_name} -PARAMETERS_NAME {sexparam} -DETECT_MINAREA {minarea} -FILTER_NAME {sexconv} -DEBLEND_NTHRESH {deb_nthresh} -DEBLEND_MINCONT {deb_mincont} -GAIN {gain} -PIXEL_SCALE {pixscale} -SEEING_FWHM {FWHM} -STARNNW_NAME {sexnnw} -BACK_TYPE AUTO -BACK_SIZE {backsize} -BACK_FILTERSIZE {backfiltersize} -BACKPHOTO_TYPE {backphoto_type} -BACKPHOTO_THICK {backphoto_thick} -CHECKIMAGE_TYPE {checkimage_type} -CHECKIMAGE_NAME {checkimage_name} -DETECT_THRESH {det_thresh}"
    if aperture is not None:
        # MAG_APER aperture diameter(s) in pixels
        # apertrue = '3.3,5.4,7.0'
        sexcom = sexcom + f" -PHOT_APERTURES {aperture}"
    print(sexcom)
    #os.system(sexcom)
    cmd    = subprocess.getstatusoutput(sexcom)
    #rms    = cmd[1].split('/')[3].split('RMS:')[-1].strip()
    #rms    = list(cmd)[-1].split('/')[4].split('RMS:')[-1].strip()

    rmsidx = str(cmd).find('RMS')
    if rmsidx == -1:
        print("SE is not running. Find the issue.")
    else :
        rms    = str(cmd)[rmsidx+5:rmsidx+12]
        secat  = ascii.read(outcat_name)
    return secat, sexcom, aperture, rms 

def calpixscale(data, hdr) :

    ny, nx = data.shape
    wcs    = WCS(hdr)

    center_x, center_y = nx / 2, ny / 2
    # 중앙 픽셀과 오른쪽 픽셀의 천체 좌표 계산
    ra_dec_center = SkyCoord.from_pixel(center_x, center_y, wcs)
    ra_dec_right = SkyCoord.from_pixel(center_x + 1, center_y, wcs)

    # 중앙 픽셀과 윗쪽 픽셀의 천체 좌표 계산
    ra_dec_top = SkyCoord.from_pixel(center_x, center_y + 1, wcs)

    # RA와 Dec에서의 각도 차이를 아크초 단위로 계산
    delta_ra = ra_dec_center.separation(ra_dec_right).arcsecond
    delta_dec = ra_dec_center.separation(ra_dec_top).arcsecond

    # 평균 픽셀 스케일 계산
    pixel_scale = (delta_ra + delta_dec) / 2  # 아크초/픽셀

    print(f"Pixel scale: {pixel_scale:.4f} arcsec/pixel")

    return pixel_scale


def process_catalog(ref, obj, ra, dec, fov, ratio):
    query_functions = {
        'PS1': {'query': iph.ps1_query, 'convert': iph.ps1_Tonry, 'suffix': 'ps1'},
        'APASS': {'query': iph.apass_query, 'convert': iph.apass_Blanton, 'suffix': 'apass'},
        'APASS10': {'query': iph.apass10_query, 'convert': iph.apass_Blanton, 'suffix': 'apass10'},
        '2MASS': {'query': iph.twomass_query, 'convert': None, 'suffix': '2mass'},
        'SMdr2': {'query': iph.SkyMapper_query, 'convert': iph.SkyMapper_Blanton, 'suffix': 'skymapper'},
        'SDSS': {'query': iph.sdss_query, 'convert': iph.sdss_Blanton, 'suffix': 'sdss'},
        'L13': {'query': iph.landolt13_query, 'convert': None, 'suffix': 'L13'},
    }
    
    if ref not in query_functions:
        raise ValueError(f"Unsupported reference catalog: {ref}")
    
    query_info = query_functions[ref]
    querycatname = f"{query_info['suffix']}-{obj}.cat"
    radius = fov / 2 * ratio if ref not in ['2MASS', 'SMdr2'] else fov * ratio
    
    if not os.path.exists(querycatname):
        print(f"Download {querycatname}...")
        if ref == 'L13':  # 특별한 옵션이 필요한 경우 처리
            querycat = query_info['query'](obj, ra, dec, radius=radius, AB=True)
        else:
            querycat = query_info['query'](obj, ra, dec, radius=radius)
        
        if query_info['convert']:
            convquerycat, convquerycat_name = query_info['convert'](querycat, obj)
        else:
            convquerycat = querycat
    else:
        print(f"{querycatname} already exists...")
        querycat = ascii.read(querycatname)
        if query_info['convert']:
            convquerycat, convquerycat_name = query_info['convert'](querycat, obj)
        else:
            convquerycat = querycat
    
    return convquerycat


def phot(imlist_name, target_catalog='', band='V', path='./lgpyred/photconf/', savecat=True, subprefix='hdCalib', ref='PS1', mykey='APER_1', sub=False, minarea=5, det_thresh=5, deb_nthresh=32, deb_mincont=0.01, backsize=128, backfiltersize = 5, backphoto_type='LOCAL', ratio=0.9, onlyzp=True, snrcut=0.1, magup=11, maglow=15):
    """
    MAG_AUTO   : AUTO 
    MAG_APER_0 : 1*FWHM
    MAG_APER_1 : 2*FWHM
    MAG_APER_2 : 3*FWHM
    MAG_APER_3 : 4*FWHM
    MAG_APER_4 : 5*FWHM
    MAG_APER_5 : 14" apeture
    MAG_APER_6 : 3"
    MAG_APER_7 : 5"
    MAG_APER_8 : 10"

    TEST
    target_catalog=''; band='V'; path='/home/somang/mylib/lgpyred/photconf/'; savecat=True; subprefix='hdCalib'; ref='PS1'; mykey='APER_1'; sub=False; minarea=5; det_thresh=3; deb_nthresh=32; deb_mincont=0.01; backsize=128; backfiltersize = 5; backphoto_type='LOCAL'; backphoto_thick=24; ratio=1; onlyzp=True; snrcut=0.1; magup=13; maglow=17.5
    """
    start_time = time.time()
    myerrkey   = 'MAGERR_'+mykey
    
    with pkg_resources.path(lgpyred.data, 'obs_spec.txt') as path:
        obscat = ascii.read( str(path))

    if savecat :
        f = open("phot.result.txt", 'w+')
        f.write(f"#{'refcat:'} {ref:<10}\n")
        f.write(f"#{'magkey:'} {mykey:<10}\n")
        #f.write(f"#{'aper_n:'} {aper_n:<10}\n")
        #f.write(f"#{'apersize:'} {str(aper_n)+'*'+mykey:<10}\n")
        f.write(f"#{'subtraction:'} {sub:<10}\n")
        f.write(f"#{'minarea:'} {minarea:<10}\n")
        f.write(f"#{'det_thresh:'} {det_thresh:<10}\n")
        f.write(f"#{'deb_nthresh:'} {deb_nthresh:<10}\n")
        f.write(f"#{'deb_mincont:'} {deb_mincont:<10}\n")
        f.write(f"#{'backsize:'} {backsize:<10}\n")
        f.write(f"#{'backfiltersize:'} {backfiltersize:<10}\n")
        #f.write(f"#{'backphoto_type:'} {backphoto_thick:<10}\n")
        f.write(f"#{'Use subimage :'} {sub}\n")
        f.write(f"{'filename'} {'obs'} {'ra'} {'dec'} {'filter'} {'date-obs'} {'mjd'} {'fwhm'} {'mag_auto'} {'magerr_auto'} {'mag'} {'magerr'} {'zp'} {'zper'} {'ul_5sig'} {'ul_3sig'} {'det'}\n") 
    imlist = glob.glob(imlist_name) ; imlist.sort()
    for inim in imlist: 
        print(f'Start photometry for {inim}')
        data, hdr   = fits.getdata(inim, header=True)
        obs_ccd     = inim.split('-')[1]
        obj         = inim.split('-')[2]
        #band        = inim.split('-')[5]
        #band        = 
        if (obs_ccd == 'BOAO_KASINICS') & (band == 'Ks'):
            band = 'K'
            
        obs_idx     = np.where(obs_ccd == obscat['obs_ccd'])[0]
        if len(obs_idx) == 0:
            raise ValueError(f"No matching obs_ccd found for {obs_ccd}")
        #binning     = float(hdr['XBINNING'])
        #pixscale    = float(obscat['pixelscale'][obs_idx].item()*binning)
        pixscale = calpixscale(data, hdr)
        gain        = float(obscat['gain'][obs_idx].item())
        fov         = obscat['fov_a'][obs_idx][0]/60.
        magl        = obscat['magl'][obs_idx]
        magh        = obscat['magh'][obs_idx]

        # Image center
        ra, dec = wcscenter(inim)
        print(f"Center coordinates: {ra}, {dec}")

        convquerycat = process_catalog(ref, obj, ra, dec, fov, ratio)
            
        # 1st SE for FWHM
        precat_name = inim[:-5]+'.pre.cat'
        print(f"==== SE config")
        print(f'pixscale : {pixscale}')
        print(f'gain     : {gain}')
        precat, presexcom, _ , _ = SE(inim, precat_name, pixscale, path=path, gain=gain, FWHM=1.2, minarea=minarea, det_thresh=det_thresh, deb_nthresh=deb_nthresh, deb_mincont=deb_mincont, backsize=backsize, backfiltersize=backfiltersize, backphoto_type=backphoto_type, backphoto_thick=24, aperture=None)

        # Pre - FWHM
        precat         = precat[(precat['MAGERR_AUTO'] < 0.2) & (precat['FLAGS'] <= np.min(precat['FLAGS'])) & (precat['MAG_AUTO'] < magl) & (precat['MAGERR_AUTO'] > magh)]
        fwhmpix        = precat['FWHM_IMAGE']
        fwhmpix        = fwhmpix[~np.isnan(fwhmpix)]

        print(fwhmpix)

        #FWHMPIX, _, FWHMPIX_std = sigma_clipped_stats(fwhmpix, sigma=5.0)
        FWHMPIX        = np.median(fwhmpix) # Gaussuan core diameter in pixel
        FWHMPIX_std    = np.std(fwhmpix)

        fwhm           = fwhmpix*pixscale

        FWHM, _, FWHM_std = sigma_clipped_stats(fwhmpix*pixscale, sigma=2.0)
        E, _, E_std       = sigma_clipped_stats(precat['ELLIPTICITY'], sigma=2.0)

        # 2nd SE
        outcat_name = inim[:-5]+'.cat'
        #aperture    = f'{aper_n*(FWHMPIX)}' # pix diameter

        magkeys = ['APER_0','APER_1', 'APER_2', 'APER_3', 'APER_4', 'APER_5', 'APER_6', 'APER_7', 'APER_8'] 
        # APER_0 - 4 : 1 - 5 * FWHM
        
        # MAG_APER aperture diameter(s) in pixels
        apertures = list(np.linspace(1, len(magkeys), num=len(magkeys)) * FWHMPIX) # in pixel diameter, -1 은 아래 유저 구경을 추가하기 위해 빼줌
        # FWHM is the diam. of the disk that contains half of the object flux.

        # Add apertures in user custom. (default.param에서 APER 개수를 임의로 추가할 것.)
        # APER_5     : Landolt aperture 14"
        # APER_6     : 3" aperture
        # APER_7     : 5" aperture
        # APER_8     : 10" aperture
        add_aper  = [14./pixscale, 3./pixscale, 5./pixscale, 10./pixscale]
        apertures = apertures + add_aper
        magkeys   = magkeys + ['APER_9', 'APER_10', 'APER_11', 'APER_12']

        aperture_dict = {key: aperture for key, aperture in zip(magkeys, apertures)}
        aperture = ','.join(str(a) for a in apertures)

        backphoto_thick = 3.*(FWHM/pixscale) # SNU IRAF photometry manual -> dannulus

        photcat, photcom, _, rms = SE(inim, outcat_name, pixscale, path=path, FWHM=FWHM, gain=gain, minarea=minarea, det_thresh=det_thresh, deb_nthresh=deb_nthresh, deb_mincont=deb_mincont, backsize=backsize, backfiltersize=backfiltersize, backphoto_type=backphoto_type, backphoto_thick=backphoto_thick, aperture=aperture)

        # Merge catalog
        merge    = iph.matching(outcat_name, convquerycat, sep = 3.)
        merge.rename_column('MAG_APER', 'MAG_APER_0')
        merge.rename_column('MAGERR_APER', 'MAGERR_APER_0')
        merge.rename_column('FLUX_APER', 'FLUX_APER_0')
        merge.rename_column('FLUXERR_APER', 'FLUXERR_APER_0')

        #if ref != 'SMdr2' :
            #MergeCoo = SkyCoord(ra=merge['ra'], dec=merge['dec'])
        #elif ref == 'SMdr2' :
            #MergeCoo = SkyCoord(ra=merge['ra']*u.deg, dec=merge['dec']*u.deg)

        try:
            MergeCoo = SkyCoord(ra=merge['ra']*u.deg, dec=merge['dec']*u.deg)
        except:
            MergeCoo = SkyCoord(ra=merge['ra'], dec=merge['dec'])

        refcat           = merge
        refcat['refnum'] = refcat['NUMBER']
        print(aperture_dict)

        ZPlist    = []
        zpstdlist = []
        lim5_list = []
        lim3_list = []

        for magkey, aper in aperture_dict.items():
            print(magkey, aper, 'pix diameter')

            if ref == 'PS1':
                #maglow = 16
                #magup = 13.5
                refcat_good      = refcat[(refcat['MAG_'+magkey] != 99) &  (refcat['MAGERR_'+magkey] < snrcut) & (refcat['FLAGS'] == 0) & (refcat[band] < maglow) & (refcat[band] > magup) ]
                
            elif ref == 'APASS' :
                #maglow = 15
                #magup = 11
                refcat_good           = refcat[(refcat['MAG_'+magkey] != 99) &  (refcat['MAGERR_'+magkey] < snrcut) & (refcat['FLAGS'] == 0) & (refcat[band] < maglow) & (refcat[band] > magup) ]
            elif ref == '2MASS' :
                if obs_ccd == 'BOAO_KASINICS':
                    # magup=11
                    # maglow = 15
                    refcat_good           = refcat[(refcat['MAG_'+magkey] != 99) &  (refcat['MAGERR_'+magkey] < snrcut) & (refcat['FLAGS'] == 0) & (refcat[band] < maglow) & (refcat[band] > magup) ]
                else :
                    # magup=11
                    # maglow = 18
                    refcat_good           = refcat[(refcat['MAG_'+magkey] != 99) &  (refcat['MAGERR_'+magkey] < snrcut) & (refcat['FLAGS'] == 0) & (refcat[band] < maglow) & (refcat[band] > magup) ]
            elif ref == 'SDSS' :
                # magup = 10
                # maglow = 16
                refcat_good           = refcat[(refcat['MAG_'+magkey] != 99) &  (refcat['MAGERR_'+magkey] < snrcut) & (refcat['FLAGS'] == 0) & (refcat[band] < maglow) & (refcat[band] > magup) ]
            elif ref == 'SMdr2' :
                # magup = 10
                # maglow = 16
                refcat_good           = refcat[(refcat['MAG_'+magkey] != 99) &  (refcat['MAGERR_'+magkey] < snrcut) & (refcat['FLAGS'] == 0) & (refcat[band] < maglow) & (refcat[band] > magup) ]
            elif ref == 'L13':
                # magup = 10
                # maglow = 15
                refcat_good           = refcat[(refcat['MAG_'+magkey] != 99) &  (refcat['MAGERR_'+magkey] < snrcut) & (refcat['FLAGS'] == 0) & (refcat[band] < maglow) & (refcat[band] > magup) ]

            if len(refcat_good) < 3:
                print('Usable star number < 3')
                refcat_good = refcat
            else:                
                pass

            # ZP 
            refcat_good['zp']   = refcat_good[band] - refcat_good[f'MAG_{magkey}']
            refcat_good['zper'] = np.sqrt(refcat_good[f'{band}err']**2 + refcat_good[f'MAGERR_{magkey}']**2) # err = refmagerr**2 (conversion error included) + SE measurement error**2

            refcat_good['refmag']    = refcat_good['zp'] + refcat_good['MAG_'+magkey]
            refcat_good['refmagerr'] = np.sqrt(refcat_good['zper']**2 + refcat_good['MAGERR_'+magkey]**2)

            #meanzp, medianzp, zpstd = sigma_clipped_stats(zp, sigma=2.0)
            if obs_ccd in ['SOAO'] :
                refcat_clip      = sigma_clip(refcat_good['zp'], sigma=2.0, maxiters=5, cenfunc=median, masked=True, copy=False)
            else :
                refcat_clip      = sigma_clip(refcat_good['zp'], sigma=2.0, maxiters=5, cenfunc=median, masked=True, copy=False)
            indx_alive       = np.where(refcat_clip.mask == False)[0]
            indx_exile       = np.where(refcat_clip.mask == True )[0]

            refcat_alive     = refcat_good[indx_alive]
            refcat_exile     = refcat_good[indx_exile]
            if magkey == mykey :
                totals, alives, exiles = refcat_good, refcat_alive, refcat_exile
            else:
                pass

            # Zero point
            ZP      = np.mean(refcat_alive['zp'])
            zpstd   = np.std(refcat_alive['zp'])
            ZPlist.append(ZP)
            zpstdlist.append(zpstd)

            # Detection limit
            limit5 = iph.limitmag(5, ZP, float(aper), float(rms))
            limit3 = iph.limitmag(3, ZP, float(aper), float(rms))
            lim5_list.append(limit5)
            lim3_list.append(limit3)

            if magkey == mykey :
                print(f"========= Photometry results {magkey} =============")
                print(inim)
                for i, m, n in zip(alives['refnum'], alives['zp'], alives['zper']) :
                    print(f"ZP of star {i} = {m:.2f}+/-{n:.2f}")
                print(f"# of alives : {len(alives)} of {len(totals)}")
                #print(f"# of masked : {len(exiles)}")
                print(f"ZP   (mean)  = {ZP:.2f}+/-{zpstd:.2f}")
                print(f"FWHM (mean)  = {FWHM:.2f}+/-{FWHM_std:.2f} arcsec")
                print(f"5sigma depth = {limit5:.2f} AB mag")
                print(f"3sigma depth = {limit3:.2f} AB mag")
                print("================================================")
                print(f"Aperture size = {aper} pix diameter")
        
        ZP_dict    = {magkey : ZPlist for ZPlist, magkey in zip(ZPlist, magkeys)}
        ZPer_dict  = {magkey : zpstdlist for zpstdlist,magkey in zip(zpstdlist, magkeys)}
        lim5_dict  = {magkey : lim5 for lim5, magkey in zip(lim5_list, magkeys)}
        lim3_dict  = {magkey : lim3 for lim3, magkey in zip(lim3_list, magkeys)}
        print('ZP is measured.')

        print(f'band     = {band}')
        print(f'ref      = {ref}')
        print(f'mykey    = {mykey}')
        print(f'myerrkey = {myerrkey}')
        #print(f'totals = {totals}')
        #print(f'alives = {alives}')
        #print(f'exiles = {exiles}')
        print(f'ZP_dict[mykey]   = {ZP_dict[mykey]}')
        print(f'ZPer_dict[mykey] = {ZPer_dict[mykey]}')
        #print(f'fwhm = {fwhm}')
        print(f'FWHM     = {FWHM}')
        print(f'FWHM_std = {FWHM_std}')
 
        PlotPhotRes(inim, band, ref, mykey, myerrkey, totals, alives, exiles, ZP_dict[mykey], ZPer_dict[mykey], fwhm, FWHM, FWHM_std)

        # Enter header
        hdr['FWHM']          = f'{FWHM:.3f}'
        hdr['FWHMSTD']       = f'{FWHM_std:.3f}'
        hdr.comments['FWHM'] = '2sigma-clipped Mean FWHM["] (SE)'
        hdr['E']             = f'{E:.3f}'
        hdr['ESTD']          = f'{E_std:.3f}'
        hdr.comments['E']    = '2sigma-clipped Mean ELLIPTICITY. 0:circle'

        hdr['mykey']            = mykey
        hdr.comments['mykey']   = 'MAG_APER key for ZP, LIM'
        hdr['PIXDIAM']          = f'{aperture_dict[mykey]:.3f}'
        hdr.comments['PIXDIAM'] = f'Aper diameter in pix ({mykey})'
        hdr['ZP']               = f'{ZP_dict[mykey]:.3f}'
        hdr['ZPER']             = f'{ZPer_dict[mykey]:.3f}'
        hdr.comments['ZP']      = 'ZP [AB], '+ref+', inmagkey='+mykey
        hdr['NTOTAL']           = f'{len(totals)}'
        hdr['NALIVE']           = f'{len(alives)}'
        hdr['NEXILE']           = f'{len(exiles)}'
        hdr.comments['NTOTAL']  = 'N of selected stars for ZP'
        hdr.comments['NALIVE']  = 'N of sigma-clipped stars for ZP'
        hdr.comments['NEXILE']  = 'N of rejected stars for ZP'
        hdr['SKYRMS']           = f'{float(rms):.3f}'
        hdr.comments['SKYRMS']  = 'skybackground rms from SE run.'
        hdr['LIM5']             = f'{lim5_dict[mykey]:.3f}'
        hdr.comments['LIM5']    = f'5sigma AB depth for a point source ({mykey})'
        hdr['LIM3']             = f'{lim3_dict[mykey]:.3f}'
        hdr.comments['LIM3']    = f'3sigma AB depth for a point source ({mykey})' 

        fits.writeto(inim, data, header=hdr, overwrite=True)

        ########### Photometry of Target (SN) ###########
        if onlyzp :
            pass
        elif onlyzp == False:
            if sub == True :
                # Photometry on Subtracted image
                part    = inim.split('-')
                part[0] = subprefix
                subim   = '-'.join(part)

                suboutcat_name = subim[:-5]+'.cat'
                photcat, photcom, _, _ = SE(subim, suboutcat_name, pixscale, gain=gain, FWHM=FWHM, minarea=minarea, det_thresh=det_thresh, deb_nthresh=deb_nthresh, deb_mincont=deb_mincont, backsize=backsize, backfiltersize=backfiltersize, backphoto_type=backphoto_type, backphoto_thick=backphoto_thick,aperture=aperture)            
                # SN target             = (154.16100585, +73.40056125)
                # Star target           = (154.3748812, +73.40496249)
            else :
                pass
            try:
                #target         = ascii.read(path+"target.cat")
                target = ascii.read(path+target_catalog)
                print(f'{path}{target_catalog}')
            except:
                print("You should decide a target to perform photometry.")
                raise Exception("No target error.")
            targetCoo      = SkyCoord(ra=target['radeg']*u.degree, dec=target['decdeg']*u.degree)

            photCoo        = SkyCoord(ra=photcat['ALPHA_J2000'], dec=photcat['DELTA_J2000'])
            target_id, target_d2d, _ = targetCoo.match_to_catalog_sky(photCoo) 

            photcat.rename_column('MAG_APER', 'MAG_APER_0')
            photcat.rename_column('MAGERR_APER', 'MAGERR_APER_0')
            photcat.rename_column('FLUX_APER', 'FLUX_APER_0')
            photcat.rename_column('FLUXERR_APER', 'FLUXERR_APER_0')

            targetcat              = photcat[target_id] 
            targetcat['sep']       = target_d2d.arcsec
            targetcat['targetnum'] = target['targetnum']
            print(targetcat)

            for i in range(len(apertures)) :
                targetcat[f'MAG_APER_{i}_ZP'], targetcat[f'MAGERR_APER_{i}_ZP'] = targetcat[f'MAG_APER_{i}'], targetcat[f'MAGERR_APER_{i}']
                targetcat[f'MAG_APER_{i}_ZP']    = ZPlist[i] +  targetcat[f'MAG_APER_{i}']
                targetcat[f'MAGERR_APER_{i}_ZP'] = np.sqrt(zpstdlist[i]**2 +  targetcat[f'MAGERR_APER_{i}']**2)

            if sub :
                targetcat.write(subim[:-5]+'.zp.cat', format='ascii', overwrite=True)
            else :
                targetcat.write(inim[:-5]+'.zp.cat', format='ascii', overwrite=True)

            for res in targetcat :
                if res['sep'] < 3.0 :
                    det_flag   = 1
                    automag    = ZP_dict[mykey] + res[f'MAG_AUTO']
                    automagerr = np.sqrt(ZPer_dict[mykey]**2 + res['MAGERR_AUTO']**2)
                    mag        = ZP_dict[mykey] + res[f'MAG_{mykey}'] 
                    magerr     = np.sqrt(ZPer_dict[mykey]**2 + res[myerrkey]**2)

                    print(f"Target {res['targetnum']}")
                    print(f"{band} = {mag:.2f}+\-{magerr:.2f} AB ({f'MAG_{mykey}'}, {obs_ccd}, {ref})")
                    print("================================================")
                    writeto = f"{inim} {obs_ccd} {targetcat['ALPHA_J2000'][0]:.7f} {targetcat['DELTA_J2000'][0]:.7f} {band} {hdr['date-obs']} {hdr['mjd']:.3f} {FWHM:.3f} {automag:.3f} {automagerr:.3f} {mag:.3f} {magerr:.3f} {ZP_dict[mykey]:.3f} {ZPer_dict[mykey]:.3f} {lim5_dict[mykey]:.3f} {lim3_dict[mykey]:.3f} {det_flag}\n"
                    print(writeto)
                    if savecat :
                        print("Write photometric result on phot.result.txt...")
                        f.write(writeto)
                    
                else :
                    print("No target is matched.")
                    det_flag   = 0
                    writeto = f"{inim} {obs_ccd} {-99:.7f} {-99:.7f} {band} {hdr['date-obs']} {hdr['mjd']:.3f} {FWHM:.3f} {-99:.3f} {-99:.3f} {-99:.3f} {-99:.3f} {ZP_dict[mykey]:.3f} {ZPer_dict[mykey]:.3f} {lim5_dict[mykey]:.3f} {lim3_dict[mykey]:.3f} {det_flag}\n"
                    print(writeto)
                    if savecat :
                        f.write(writeto)      
    if savecat :
        f.close()
        #os.system('pluma '+"phot.result.txt &")
    end_time = time.time()
    print(f'*** run time : {end_time - start_time:.2f}s')
    if onlyzp :
        return None
        #return merge
    elif onlyzp == False:
        return targetcat


def subphot(imlist_name='Cal*com.fits', target_catalog='SN2025dr.cat', band='i', path='/home/lim9/miniconda3/lib/python3.9/site-packages/lgpyred/photconf/', mykey='APER_1', subprefix='hdCalib', minarea=3, det_thresh=3, deb_nthresh=32, deb_mincont=0.01,  backsize=128, backfiltersize=5,  backphoto_type='LOCAL', savecat=True):
    """
    ZP가 헤더에 존재하는 이미지에서 특정 구경(APER_1)으로 측광을 수행하는 함수

    Parameters:
    imlist_name (str): 이미지 파일 패턴 (예: 'Cal*.fits')
    target_catalog (str): 측광 대상(예: 초신성)의 좌표가 담긴 파일
    path (str): SExtractor 설정 파일 경로
    mykey (str): 사용할 구경 (디폴트: 'APER_1')
    subprefix (str): 보정된 이미지의 접두사
    minarea, det_thresh, deb_nthresh, deb_mincont, backsize, backfiltersize, backphoto_type: SE 설정

    Returns:
    Table: 측광 결과 테이블
    """
    start_time = time.time()

    imlist = glob.glob(imlist_name); imlist.sort()

    if not imlist:
        print("No matching images found.")
        return None

    try:
        target = ascii.read(path + target_catalog)
    except FileNotFoundError:
        print(f"Target catalog {target_catalog} not found. Exiting.")
        return None

    with open('phot_result.txt', "+w") as f:
        f.write("filename obs ra dec filter date-obs mjd fwhm mag_auto magerr_auto mag magerr zp zper ul_5sig ul_3sig det\n")

        for inim in imlist:
            print(f'Processing {inim}...')
            obs_ccd     = inim.split('-')[1]
            obj         = inim.split('-')[2]

            # FITS 헤더에서 ZP 읽기
            data, hdr = fits.getdata(inim, header=True)
            try:
                ZP    = float(hdr['ZP'])
                e_ZP = float(hdr['ZPER'])
                FWHM  = float(hdr['FWHM'])
                lim3  = float(hdr['lim3'])
                lim5  = float(hdr['lim5'])
            except KeyError:
                print(f"ZP 정보가 {inim}에 없음. 건너뜀.")
                continue

            # Pixel Scale 계산
            pixscale = calpixscale(data, hdr)

            # Gain 정보
            gain = hdr.get('EGAIN', 1)

            # Photometry

            subim = inim.replace("Calib", subprefix)

            # Aperture 설정
            apertures = np.array([(1.*FWHM)/pixscale, (2.*FWHM)/pixscale, (3.*FWHM)/pixscale, (4.*FWHM)/pixscale, (5.*FWHM)/pixscale, 14./pixscale, 3./pixscale, 5./pixscale, 10./pixscale])
            
            aperture = ",".join(apertures.astype(str))
            outcat_name = subim[:-5] + '.zp.cat'

            backphoto_thick = 3. * (FWHM/pixscale)
            photcat, _, _, _ = SE(
                subim, outcat_name, pixscale, path=path, gain=gain, FWHM=FWHM,
                minarea=minarea, det_thresh=det_thresh, deb_nthresh=deb_nthresh, deb_mincont=deb_mincont,
                backsize=backsize, backfiltersize=backfiltersize, backphoto_type=backphoto_type,
                backphoto_thick=backphoto_thick, aperture=aperture
            )

            photcat.rename_column('MAG_APER', 'MAG_APER_0')
            photcat.rename_column('MAGERR_APER', 'MAGERR_APER_0')
            photcat.rename_column('FLUX_APER', 'FLUX_APER_0')
            photcat.rename_column('FLUXERR_APER', 'FLUXERR_APER_0')

            # APER_1 필드 확인
            if f'MAG_{mykey}' not in photcat.colnames or f'MAGERR_{mykey}' not in photcat.colnames:
                print(f"{mykey} 데이터가 {inim}에 없음. 건너뜀.")
                continue
            # 목표 천체(예: 초신성) 좌표 불러오기
            targetCoo = SkyCoord(ra=target['radeg'] * u.degree, dec=target['decdeg'] * u.degree)
            photCoo = SkyCoord(ra=photcat['ALPHA_J2000'], dec=photcat['DELTA_J2000'])
            target_id, target_d2d, _ = targetCoo.match_to_catalog_sky(photCoo)

            # 가장 가까운 소스 선택
            targetcat = photcat[target_id]
            targetcat['sep'] = target_d2d.arcsec
            targetcat['targetnum'] = target['targetnum']

            # ZP 적용하여 측광값 변환
            targetcat[f'MAG_{mykey}_ZP'] = ZP + targetcat[f'MAG_{mykey}']
            targetcat[f'MAGERR_{mykey}_ZP'] = np.sqrt(e_ZP**2 + targetcat[f'MAGERR_{mykey}']**2)

            # 저장할 파일명 결정
            #targetcat.write(output_filename, format='ascii', overwrite=True)

            # 결과 파일에 저장
            for res in targetcat :
                if res['sep'] < 3.0 :
                    det_flag   = 1
                    automag    = ZP + res[f'MAG_AUTO']
                    automagerr = np.sqrt(e_ZP**2 + res['MAGERR_AUTO']**2)
                    mag        = ZP + res[f'MAG_{mykey}'] 
                    magerr     = np.sqrt(e_ZP**2 + res[f'MAGERR_{mykey}']**2)

                    print(f"Target {res['targetnum']}")
                    print(f"{band} = {mag:.2f}+\-{magerr:.2f} AB ({f'MAG_{mykey}'}, {obs_ccd})")
                    print("================================================")
                    writeto = f"{inim} {obs_ccd} {targetcat['ALPHA_J2000'][0]:.7f} {targetcat['DELTA_J2000'][0]:.7f} {band} {hdr['date-obs']} {hdr['mjd']:.3f} {FWHM:.3f} {automag:.3f} {automagerr:.3f} {mag:.3f} {magerr:.3f} {ZP:.3f} {e_ZP:.3f} {lim5:.3f} {lim3:.3f} {det_flag}\n"
                    print(writeto)
                    if savecat :
                        print("Write photometric result on phot.result.txt...")
                        f.write(writeto)
                    
                else :
                    print("No target is matched.")
                    det_flag   = 0
                    writeto = f"{inim} {obs_ccd} {-99:.7f} {-99:.7f} {band} {hdr['date-obs']} {hdr['mjd']:.3f} {FWHM:.3f} {-99:.3f} {-99:.3f} {-99:.3f} {-99:.3f} {ZP:.3f} {e_ZP:.3f} {lim5:.3f} {lim3:.3f} {det_flag}\n"
                    print(writeto)
                    if savecat :
                        f.write(writeto)   

    end_time = time.time()
    print(f'*** Run time : {end_time - start_time:.2f}s')

    print(f"Photometry results saved to phot_result.txt")




