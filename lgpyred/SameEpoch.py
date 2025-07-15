def SameEpoch(imlist_name, sep=5.):
    """
    Return a list of images taken at the same epoch (< Separation).
    Separation default is 5 minutes.
    """
    import glob
    import numpy as np
    from astropy.io import fits
    from astropy.table import Table

    # Gather image files
    imlist = glob.glob(imlist_name)
    imlist.sort()

    if not imlist:
        print('No images found.')
        return []

    print('Separation is : {} Minutes.'.format(sep))

    # Initialize lists to store information
    image_data = []

    # Extract header information from each image
    for inim in imlist:
        hdr = fits.getheader(inim)
        mjd = hdr['MJD']
        exptime = hdr['EXPTIME']
        objectname = inim.split('-')[2]
        band = inim.split('-')[5]
        image_data.append((inim, objectname, band, mjd, exptime))

    # Sort by MJD
    image_data.sort(key=lambda x: x[3])

    # Initialize variables for grouping
    groups = []
    current_group = [image_data[0][0]]
    last_objectname = image_data[0][1]
    last_band = image_data[0][2]
    last_mjd = image_data[0][3]

    # Group images by objectname, band, and separation time
    for i in range(1, len(image_data)):
        inim, objectname, band, mjd, exptime = image_data[i]
        
        # Check if the current image should start a new group
        if objectname != last_objectname or band != last_band or (mjd - last_mjd) * 24 * 60 > sep:
            groups.append(current_group)
            current_group = [inim]
        else:
            current_group.append(inim)
        
        # Update last seen values
        last_objectname = objectname
        last_band = band
        last_mjd = mjd
    
    # Append the last group
    if current_group:
        groups.append(current_group)

    return groups


'''
def SameEpoch(imlist_name, sep=5.):
	"""
	Return a list of images taken at the same epoch (< Separation).
	Separation default is 5 minutes.
	"""
	import glob
	import numpy as np 
	from astropy.io import fits
	from astropy.table import Table
	imlist = glob.glob(imlist_name); imlist.sort()
	print('Separation is : '+str(sep)+' Minutes.')
	camera	     = imlist[0].split('-')[1]
	object_name  = imlist[0].split('-')[2]
	band		 = imlist[0].split('-')[5]
	exptime = []
	image   = []
	mjd  	= []
	for i in range(len(imlist)) :
		inim = imlist[i]
		hdr  = fits.getheader(inim)
		mjd.append(hdr['MJD'])
		exptime.append(hdr['exptime'])
		image.append(inim)
	tbl = Table({'image' : image, 'mjd': mjd, 'exptime':exptime}, names = ['image', 'mjd', 'exptime'])
	epoch_prev = []
	for j in range(len(tbl['mjd'])):
		dmjd_dum = np.array(tbl['mjd']) - tbl['mjd'][j]
		idx = np.where(abs(dmjd_dum*24*60) < sep)[0]
		epoch_prev.append(list(tbl['image'][idx]))
	res = list(set(tuple(sorted(sub)) for sub in epoch_prev))
	res.sort()
	return res
'''
def SameEpoch2(imlist_name, sep=5.):
	"""
	Return a list of images taken at the same epoch (< Separation).
	Separation default is 5 minutes.
	"""
	import glob
	import numpy as np 
	from astropy.io import fits
	from astropy.table import Table
	imlist = glob.glob(imlist_name)
	imlist.sort()
	print('Separation is : '+str(sep)+' Minutes.')
	camera	     = imlist[0].split('-')[1]
	object_name  = imlist[0].split('-')[2]
	band		 = imlist[0].split('-')[5]
	exptime = []
	image   = []
	mjd  	= []
	for i in range(len(imlist)) :
		inim = imlist[i]
		hdr  = fits.getheader(inim)
		mjd.append(hdr['MJD'])
		exptime.append(hdr['exptime'])
		image.append(inim)
	tbl = Table({'image' : image, 'mjd': mjd, 'exptime':exptime}, names = ['image', 'mjd', 'exptime'])
	epoch_prev = []
	for j in range(len(tbl['mjd'])):
		dmjd_dum = np.array(tbl['mjd']) - tbl['mjd'][j]
		idx = np.where(abs(dmjd_dum*24*60) < sep)[0]
		epoch_prev.append(list(tbl['image'][idx]))
	res = list(set(tuple(sorted(sub)) for sub in epoch_prev))
	res.sort()
	res2 = []
	for j in range(len(res)) :
		res2.append(list(res[j]))
	return res2

def ImageGroup(imlist_list='Remap*0.fits', N=5):
	import glob
    # imlist에 해당하는 파일 이름들을 패턴에 맞게 필터링
	imlist = glob.glob(imlist_list); imlist.sort()
    
    # 이미지 개수가 홀수인 경우 마지막 리스트에 남는 이미지 추가
	if len(imlist) % N != 0:
		imlist.append(None)
    
	#인접한 이미지를 N개씩 묶어주기
	grouped_images = [imlist[i:i + N] for i in range(0, len(imlist), N)]
    
	return grouped_images
