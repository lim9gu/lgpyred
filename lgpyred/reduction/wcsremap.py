def wcsremap(srcim, tmpim, outim, verbose=False):
	'''
	1. Description
	run wcsremap-1.0.1 by Andrew Becker 
	wcsremap -template template.fits -source input.fits -outIm input_remapped.fits 

	srcim	   (Source)   : Larger  pixscale -> no change
	tmpim	   (Template) : Smaller pixscale -> change

	Small pixscale image ==> Large pixscale image
	'''
	import os
	import glob
	if verbose == True :
		com		= 'wcsremap -template '+tmpim+' -source '+srcim+' -outIm '+outim+' -v'
	elif verbose == False :
		com		= 'wcsremap -template '+tmpim+' -source '+srcim+' -outIm '+outim
	os.system(com)
	print(outim)

def run_wcsremap(imlist_name, refim, large=True):
	"""
	large = True : Sci pix is larger or the same with ref.
	large = False : Sci pix is smaller or when you make stacking image of one epoch. refim should be one of the image set, not other reference.
	"""
	import os
	import glob
	imlist  = glob.glob(imlist_name); imlist.sort()
	#obs_sci   = imlist[0].split('-')[1]
	#obs_ref   = refim.split('-')[1]
	if large == False :
		# pixscale of sciim is small than refim
		# sciim will be controlled.
		print('Pixscale_original < Pixscale_reference')
		for i in range(len(imlist)):
			inim = imlist[i]
			parts	= inim.split('-')
			parts[0]= 'Remap'
			outim	= '-'.join(parts)	
			os.system('rm '+outim)
			wcsremap(inim, refim, outim, verbose=False)
	elif large == True :
		# pixscale of sciim is larger than refim
		# refim will be controlled.
		print('Pixscale_original >= Pixscale_reference')
		for i in range(len(imlist)):
			inim = imlist[i]
			parts	= inim.split('-')
			parts[0]= 'Remap'
			outim	= '-'.join(parts)
			os.system('rm '+outim)
			wcsremap(refim ,inim, outim, verbose=False)
	print('Done.')

def run_wcsremap_wrap(imlist_name, refim, large=True, wcsregister=True):
	"""
	Use wcsremap wrapper
	main(tmp, src, out, wcsregister=True)
	tmp가 src로 맞춰져서 out으로 나옴.
	large = True : Sci pix is larger or the same with ref.
	large = False : Sci pix is smaller or when you make stacking image of one epoch. refim should be one of the image set, not other reference.
	"""
	import os
	import glob
	from lgpyred.reduction.wcsremap_wrap import main as remap
	imlist  = glob.glob(imlist_name); imlist.sort()
	if large == False :
		# pixscale of sciim is small than refim
		# sciim will be controlled.
		print('Pixscale_original < Pixscale_reference')
		for i in range(len(imlist)):
			inim = imlist[i]
			parts	= inim.split('-')
			parts[0]= 'Remap'
			outim	= '-'.join(parts)	
			os.system('rm '+outim)
			remap(inim, refim, outim, wcsregister=wcsregister)
	elif large == True :
		# pixscale of sciim is larger than refim
		# refim will be controlled.
		print('Pixscale_original >= Pixscale_reference')
		for i in range(len(imlist)):
			inim = imlist[i]
			parts	= inim.split('-')
			parts[0]= 'Remap'
			outim	= '-'.join(parts)
			os.system('rm '+outim)
			remap(refim ,inim, outim, wcsregister=wcsregister)
	print('Done.')
	return outim

def run_wcsremap_epoch(imlist_name, sep=5., wcsregister=True):
	"""
	large = True : Sci pix is larger or the same with ref.
	large = False : Sci pix is smaller or when you make stacking image of one epoch.
	"""
	import os, glob
	from lgpyred.SameEpoch import SameEpoch
	from lgpyred.reduction.wcsremap import run_wcsremap_wrap
	res = SameEpoch(imlist_name, sep=sep)
	for i in range(len(res)):
		imlist = list(res[i])
		for j in range(len(imlist)):
			inim = imlist[j]
			refim = imlist[0]
			run_wcsremap_wrap(inim, refim, large=False, wcsregister=wcsregister)
	print('Done.')
