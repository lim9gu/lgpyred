from astropy.io import ascii
import importlib.resources as pkg_resources
import lgpyred.data  

def wcscenter(inim) :
	"""
	Measure image WCS center coordinates.
	"""
	from astropy.wcs import WCS
	from astropy.io import fits
	data, hdr = fits.getdata(inim, header=True)
	w		 = WCS(hdr)
	NAXIS1, NAXIS2 = hdr['NAXIS1'], hdr['NAXIS2']
	racen, deccen = w.all_pix2world(NAXIS1/2.,NAXIS2/2.,0)
	return racen.base[0][0], deccen.base[0][1]

def hdrcheck(imlist_name, ccd, fov, racen='', deccen=''):
	"""
	Check header.
	Put galaxy name comparing to IMSNG catalog. 
	Put JD, MJD information into the image header.

	SN2019ein 208.4484215 40.29569905

	"""
	import os
	import sys
	import glob
	import astropy.units as u
	from astropy.io import fits
	from astropy.io import ascii
	from astropy.time import Time	
	from astropy.coordinates import SkyCoord

	with pkg_resources.path(lgpyred.data, 'alltarget.dat') as path:
		all_catname = str(path)
		allcat = ascii.read(all_catname)

	ra, dec = all_cat['ra'], all_cat['dec']
	radeg, decdeg = [], []
	i=0
	for i in range(len(all_cat)) :
		c   = SkyCoord(str(ra[i])+' '+str(dec[i]), unit=(u.hourangle, u.deg))
		radeg.append(c.ra.deg)
		decdeg.append(c.dec.deg)
	all_cat['radeg'] = radeg
	all_cat['decdeg'] = decdeg
	coo_all = SkyCoord(radeg, decdeg, unit=(u.deg,u.deg))

	imlist = glob.glob(imlist_name)
	imlist.sort()
	i=0
	for i in range(len(imlist)) :
		inim			= imlist[i]
		print(inim)
		data, hdr	   = fits.getdata(inim, header=True)
		mjd0   = 2400000.5
		if ccd == 'MAO_FLI' :
			racen,deccen = wcscenter(inim)
			t		  = Time(hdr['DATE-OBS'], format='isot', scale='utc')
			jd		  = t.jd
			mjd	      = t.mjd
			hdr['JD'] = round(jd, 5)
		if ccd == 'MAO_SNUCAM' :
			racen,deccen = wcscenter(inim)
			t		  = Time(hdr['utdate']+'T'+hdr['utstart'], format='isot', scale='utc')
			jd		  = t.jd
			mjd	      = t.mjd
			hdr['JD'] = round(jd, 5)
		elif ccd == 'KHAO_MDFTS' :
			racen, deccen = wcscenter(inim)
			t         = Time(hdr['DATE-OBS'], format='isot', scale='utc')
			jd		  = t.jd
			mjd	      = t.mjd
		elif ccd == 'MCD30INCH':
			racen,deccen = wcscenter(inim)
			print(hdr['DATE'])
			t		  = Time(hdr['DATE'], format='isot', scale='utc')
			print(t)
			jd		  = t.jd
			mjd	      = t.mjd
			hdr['JD'] = round(jd, 5)
			hdr['filter'] = inim.split('-')[5]	
		elif ccd == 'LOAO_E2V' :
			racen,deccen = wcscenter(inim)
			t		  = Time(hdr['DATE'], format='isot', scale='utc')
			jd		  = t.jd
			mjd	      = t.mjd
			hdr['JD'] = round(jd, 5)
		elif ccd == 'UKIRT' :
			#if (racen == '') & (deccen == '') :
			#	print('Please specify image center (ra, dec) keyword.')
			band          =  inim.split('-')[5]
			hdr['filter'] = band
			#racen   = float(racen)
			#deccen  = float(deccen)
			racen,deccen = wcscenter(inim)
			utdate  = inim.split('-')[3]
			utdate  = utdate[0:4]+'-'+utdate[4:6]+'-'+utdate[6:8]
			utstart = inim.split('-')[4]
			utstart = utstart[0:2]+':'+utstart[2:4]+':'+utstart[4:6]
			t  = Time(utdate+'T'+utstart, format='isot', scale='utc')
			jd = t.jd
			hdr['date-obs'] = utdate+'T'+utstart
			hdr['JD'] = jd
			mjd = jd - mjd0
		elif ccd != 'MAO_SNUCAM' : 
			racen, deccen = wcscenter(inim)
			jd   = hdr['JD']
			mjd  = jd - mjd0
		hdr['MJD']	  = round(mjd,5)
		coo_target	  = SkyCoord(racen, deccen, unit=(u.deg, u.deg))
		indx, d2d, d3d  = coo_target.match_to_catalog_sky(coo_all)
		if d2d.arcmin > fov/2. :
			print('Coordinates of the image are not in IMSNG catalog. No matching. Maybe you obtained wrong field. OR Non-IMSNG target.' )
			fits.writeto(inim, data, header=hdr, overwrite=True)
			print('Only MJD is entered in image header.')
			pass
		elif d2d.arcmin < fov/2. :
			obj  = all_cat[indx]['obj']
			print('======================================')
			print(obj+ ' is matched.')
			print(str(round(d2d.arcmin[0],3))+ ' arcmin apart')
			print('======================================')
			hdr['object']   = obj
			fits.writeto(inim, data, header=hdr, overwrite=True)
	print('Header info inspection is finished.')
