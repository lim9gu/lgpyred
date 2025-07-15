# Image combine series
import glob
import os, sys
import subprocess
import numpy as np 
from pyraf import iraf
from astropy.io import fits
from astropy.time import Time
from lgpyred.SameEpoch import SameEpoch

def imcombine_set(imlist_name, combine='median', reject='sigclip'):
    """
    Combine specific image set and make stacked image.
    > imcombine_set('DfCal*20190502-05*.fits')
    MJD will be added using median value of image set.

    none      - No rejection
    minmax    - Reject the nlow and nhigh pixels
    ccdclip   - Reject pixels using CCD noise parameters
    crreject  - Reject only positive pixels using CCD noise parameters
    sigclip   - Reject pixels using a sigma clipping algorithm
    avsigclip - Reject pixels using an averaged sigma clipping algorithm
    pclip     - Reject pixels using sigma based on percentiles
    """
    imlist = glob.glob(imlist_name); imlist.sort()
    #mjd, exptime, date = [], [], []
    camera = imlist[0].split('-')[1]

    Exp_list  = []
    MJD_list  = []
    DATEOBS_list = []
    UT_list = []

    for im in imlist :
        hdr = fits.getheader(im)
        if ('mjd' in hdr) | ('MJD' in hdr) :
            try:
                MJD_list.append(hdr['MJD'])
            except:
                MJD_list.append(hdr['mjd'])
        elif ('MJD' not in hdr) & ('DATE-OBS' in hdr ) :
            DATEOBS_list.append(hdr['DATE-OBS'])
        elif ('MJD' not in hdr) & ('UTDATE' in hdr) :
            UT_list.append(hdr['UTDATE']+'T'+hdr['UTSTART']) 
        else:
            print('No observation time.')
        Exp_list.append(float(hdr['exptime']))
    
    NCOMBINE = len(imlist)
    INTEXP   = np.sum(Exp_list)
    if MJD_list != [] :
        MEAN_MJD = np.mean(MJD_list)
    elif DATEOBS_list != [] :
        MJD_list = []
        for date in DATEOBS_list :
            t = Time(date, scale = 'utc', format='isot')
            MJD_list.append(t.mjd)
        MEAN_MJD = np.mean(MJD_list)
    elif UT_list != [] :
        MJD_list = []
        for ut in UT_list :
            t = Time(date, scale = 'utc', format='isot')
            MJD_list.append(t.mjd)
        MEAN_MJD = np.mean(MJD_list)        
    COM_MJD  = MEAN_MJD
    COM_time = Time(COM_MJD, scale = 'utc', format='mjd')
    COM_DATE = COM_time.fits[0:-4]
    part0 = imlist[0].split('-')[0]
    print('MJD of combined image = '+str(COM_MJD))
    print('UT of combined image = '+str(COM_DATE))
    os.system('ls '+imlist_name+' > comb.list')
    object_name = hdr['object'] 
    if camera == 'LOAO_E2V' :
        band = hdr['filter'][0]
    else:
        band = hdr['filter']
    newhdr = hdr
    newhdr['MJD']      = COM_MJD
    newhdr['EXPTIME']  = INTEXP 
    newhdr['DATE-OBS'] = COM_DATE
    newhdr['NCOMBINE'] = NCOMBINE
    newhdr['UTDATE']   = COM_DATE.split('T')[0]
    newhdr['UTSTART']  = COM_DATE.split('T')[1]

    #newimage = 'com'+part0+'-'+camera+'-'+object_name+'-'+str(newhdr['UTDATE'][0:4])+str(newhdr['UTDATE'][5:7])+str(newhdr['UTDATE'][8:10])+'-'+str(newhdr['UTSTART'][0:2])+str(newhdr['UTSTART'][3:5])+str(newhdr['UTSTART'][6:8])+'-'+band+'-'+str(int(INTEXP))+'.fits'
    #newimage = f"{part0}-{camera}-{object_name}-{str(newhdr['UTDATE'][0:4])}{newhdr['UTDATE'][5:7]}{str(newhdr['UTDATE'][8:10])}-{str(newhdr['UTSTART'][0:2])}{str(newhdr['UTSTART'][3:5])}{str(newhdr['UTSTART'][6:8])}-{band}-{str(int(INTEXP))}.com.fits"
    newimage = f"Calib-{camera}-{object_name}-{str(newhdr['UTDATE'][0:4])}{newhdr['UTDATE'][5:7]}{str(newhdr['UTDATE'][8:10])}-{str(newhdr['UTSTART'][0:2])}{str(newhdr['UTSTART'][3:5])}{str(newhdr['UTSTART'][6:8])}-{band}-{str(int(INTEXP))}.com.fits"
    print(newimage)
    iraf.imcombine('@comb.list', newimage, combine = combine, reject=reject, scale='none', zero='mode')
    newdata = fits.getdata(newimage)
    fits.writeto(newimage, newdata, header=newhdr, overwrite=True)
    print('New image and header is entered.')


def imcombine_epoch(imlist_name, sep=5., combine='median', reject='sigclip'):
    """
    Make images taken in the same epoch into one set and combine.
    You can control the separation. Default is 5 minutes.

    Parameters:
    imlist_name : str
        The pattern to match FITS files.
    sep : float
        Time separation in minutes to group images into the same epoch.
    reject : str
        The rejection algorithm to use during combination.
    """
    # Group images by epoch
    res = SameEpoch(imlist_name, sep=sep)

    # Process each group of images
    for image_group in res:
        camera      = image_group[0].split('-')[1]
        object_name = image_group[0].split('-')[2]
        band        = image_group[0].split('-')[5]

        print(image_group)
        print('will be set as the same epoch.')

        # Extract MJD and exposure time from each image
        mjds = [float(fits.getheader(img)['mjd']) for img in image_group]
        exptimes = [float(fits.getheader(img)['exptime']) for img in image_group]

        # Compute new combined properties
        total_exp   = np.sum(exptimes)
        mean_mjd    = np.mean(mjds)
        t           = Time(mean_mjd, format='mjd')
        new_dateobs = t.isot
        print("MJD values: ", mjds)
        print("Mean MJD: ", mean_mjd)
        print("New date obs (ISO): ", new_dateobs)
        yy, mm, dd  = t.isot[:4], t.isot[5:7], t.isot[8:10]
        hh, min, ss = t.isot[11:13], t.isot[14:16], t.isot[17:19]

        new_image_name = f'Calib-{camera}-{object_name}-{yy}{mm}{dd}-{hh}{min}{ss}-{band}-{int(total_exp)}.com.fits'

        # Combine images using iraf.imcombine
        iraf.imcombine(','.join(image_group), new_image_name, combine=combine, reject=reject, scale='none', zero='mode')


'''
def imcombine_epoch(imlist_name, sep=5., reject='sigclip'):
	"""
	Make images taken in the same epoch into one set and combine.
	You can control the separation. Default is 5 minutes.

	imlist_name = 'Remap*.fits'
	part0 = Calib, Remap,... front name of combined image
	sep = 5. minutes
	"""
	import glob
	import os, sys
	import numpy as np 
	from pyraf import iraf
	from lgpyred.SameEpoch import SameEpoch
	from astropy.io import fits
	from astropy.time import Time
	from astropy.table import Table

	res = SameEpoch(imlist_name, sep=sep)
	for ii in range(len(res)) :
		image_to_com = list(res[ii])
		camera       = image_to_com[0].split('-')[1]
		object_name  = image_to_com[0].split('-')[2]
		band         = image_to_com[0].split('-')[5]
		print(image_to_com)
		print('will be set as the same epoch.')
		mjd_to_com   = []
		exp_to_com   = []
		for iii in range(len(image_to_com)) :
			mjd_to_com.append(float(fits.getheader(image_to_com[iii])['mjd']))
			exp_to_com.append(float(fits.getheader(image_to_com[iii])['exptime']))
		new_exp     = np.sum(exp_to_com)
		mean_mjd	= np.mean(mjd_to_com)
		t			= Time(mean_mjd, format = 'mjd')
		new_dateobs = t.isot
		yy, mm, dd  = t.isot[0:4], t.isot[5:7], t.isot[8:10]
		hh, min, ss = t.isot[11:13], t.isot[14:16],  t.isot[17:19]

		newimage	 = 'Calib-'+camera+'-'+object_name+'-'+str(yy)+str(mm)+str(dd)+'-'+str(hh)+str(min)+str(ss)+'-'+band+'-'+str(int(new_exp))+'.com.fits'
		iraf.imcombine(','.join(image_to_com), newimage, combine = 'median', reject=reject, scale='none', zero='mode')
'''
def imcombine_group(imlist_name, N=5, reject='sigclip'):
	"""
	Make images taken in the same epoch into one set and combine.
	You can control the separation. Default is 5 minutes.

	imlist_name = 'Remap*.fits'
	part0 = Calib, Remap,... front name of combined image
	N = 5. images
	"""
	import glob
	import os, sys
	import numpy as np 
	from pyraf import iraf
	from lgpyred.SameEpoch import ImageGroup
	from astropy.io import fits
	from astropy.time import Time
	from astropy.table import Table

	res = ImageGroup(imlist_name, N=N)
	for ii in range(len(res)) :
		image_to_com = list(res[ii])
		camera = image_to_com[0].split('-')[1]
		object_name = image_to_com[0].split('-')[2]
		band = image_to_com[0].split('-')[5]
		print(image_to_com)
		print('will be set as the same epoch.')
		mjd_to_com   = []
		exp_to_com   = []
		for iii in range(len(image_to_com)) :
			mjd_to_com.append(float(fits.getheader(image_to_com[iii])['mjd']))
			exp_to_com.append(float(fits.getheader(image_to_com[iii])['exptime']))
		new_exp     = np.sum(exp_to_com)
		mean_mjd	= np.mean(mjd_to_com)
		t			= Time(mean_mjd, format = 'mjd')
		new_dateobs = t.isot
		yy, mm, dd  = t.isot[0:4], t.isot[5:7], t.isot[8:10]
		hh, min, ss = t.isot[11:13], t.isot[14:16],  t.isot[17:19]

		newimage	 = 'Calib-'+camera+'-'+object_name+'-'+str(yy)+str(mm)+str(dd)+'-'+str(hh)+str(min)+str(ss)+'-'+band+'-'+str(int(new_exp))+'.com.fits'
		iraf.imcombine(','.join(image_to_com), newimage, combine = 'median', reject=reject, scale='none', zero='mode')