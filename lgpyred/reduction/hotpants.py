import glob, os
import numpy as np
from astropy.io import fits, ascii
import importlib.resources as pkg_resources
import lgpyred.data
'''
def hotpants(sciim, refim, ngflag=False, sigma_match='', conv='t', scale='i', iu=5000000, tu=5000000, il=-10000, tl=-10000, igain=1.33, tgain=1.26, irdnoise=12.6, trdnoise=361., ssf='') :
	"""
    hotpants running in TARS
	"""

	outfile  ='hd'+sciim
	convfile ='hc'+sciim

	# setup

	# (1) Stamp, substamp

	nsx   = 12 # 이미지의 X축과 Y축에서 스탬프의 개수를 설정합니다. 기본값은 각각 10입니다.스탬프는 이미지의 특정 구역으로 나누어진 작은 영역이며, 각 영역에서 PSF와 커널을 개별적으로 계산합니다. (10) 
	nsy   = 12 # (10)

	nss   = 3 # 각 스탬프 내에서 사용할 서브스탬프 수를 설정합니다. 적절한 수를 설정해 문제 영역을 최소화할 수 있습니다. (3) 
	rss   = 15 # 각 스탬프의 반지름을 조정하여, 더 작은 영역을 검사하고 문제가 있는 부분을 피할 수 있습니다. (15) 


	# (2) Sigma clipping
	
	ks    = 2 # 커널 피팅 과정에서 이상치를 제거하기 위해 사용됩니다. 기본적으로 2.0으로 설정되며, 이 값은 커널 내 픽셀 값이 평균으로부터 몇 시그마 이상 벗어날 경우 이를 이상치로 간주하고 제거하는 기준입니다. (2)

	ssig  = 3 # -ks는 커널 피팅에서, -ssig는 통계적 분석에서 Sigma 클리핑 임계값을 설정하는 옵션입니다. 이 값을 통해 이상치나 노이즈가 있는 스탬프를 걸러낼 수 있으며, 이를 적절히 설정하지 않으면 NaN값이 발생할 수 있습니다. 이미지의 노이즈 특성에 맞게 최적의 값을 찾아야 합니다. (3.0) 

	# (3) Kernel order & Composition
	
	#r     = 10 # 컨볼루션 커널의 반지름을 설정하는 옵션입니다.커널의 크기는 이미지의 PSF 크기에 맞춰 조정해야 합니다. 너무 작은 값을 설정하면 PSF의 전체 범위를 포착하지 못하고, 너무 큰 값을 설정하면 불필요한 영역이 포함되어 노이즈가 증가할 수 있습니다.(10)

	ko    = 2 # 커널이 이미지 내에서 공간적으로 어떻게 변화하는지를 설명하는 다항식의 차수입니다.PSF가 이미지 내에서 변화할 때, 이 옵션을 적절히 설정하여 복잡한 PSF 변화를 반영할 수 있습니다. 잘못된 차수 설정은 PSF를 부정확하게 모델링할 수 있으므로, 상황에 맞게 조정하는 것이 중요합니다.(2) 

	bgo   = 1 # 이미지 내에서 배경광의 변화를 설명하는 다항식의 차수입니다. 배경이 균일하지 않은 이미지에서 적절한 값을 설정하면 배경 차이를 효과적으로 보정할 수 있습니다. 배경 변화를 정확하게 반영하지 않으면 차이 이미지에 잔여 패턴이나 거짓 신호가 남을 수 있습니다. (1) 

	kfm   = 0.99 #커널의 절대 값 합 중 유효한 값을 가질 수 있는 픽셀의 비율을 설정하는 옵션이에요. 기본값은 0.990로, 전체 커널의 99%를 유효한 값으로 간주합니다.

	# ng 

	# (4) RMS Threshold and Fit Threshold

	ft    =  21 # 커널 피팅에서 유효한 중심값을 결정하기 위한 RMS 임계값을 설정합니다. 기본값은 20으로 설정되어 있습니다. 이 값은 스탬프의 중심값(centroid)이 유효한지 여부를 결정하며, 특정 스탬프에서 이 임계값을 넘으면 해당 스탬프가 이상치로 간주될 수 있습니다. (20.0)	

	# ft = 7 # DOAO 필드 안에 별이 X, MAO_FLI도 별이 X
	afssc = 1 # 스탬프의 자동 중심 탐색을 활성화하거나 비활성화합니다.PSF가 잘 정의된 영역을 찾고, 스탬프의 중심을 자동으로 선택하여 PSF가 왜곡된 영역을 피할 수 있도록 설정하는 것이 중요합니다. (1)

	#omi   = 'mk'+sciim

	#savexy = 'st'+sciim

	if ngflag == False :
		ng = '3 6 0.70 4 1.50 2 3.00'
		
		com = f'hotpants -c {conv} -n {scale} -iu {iu} -iuk {iu} -tu {tu} -tuk {tu} -il {il} -tl {tl} -v 0 -inim {sciim} -tmplim {refim} -outim {outfile} -oci {convfile} -ks {ks} -ko {ko} -bgo {bgo} -ft {ft} -ssig {ssig} -ig {igain} -tg {tgain} -ir {irdnoise} -tr {trdnoise} -afssc {afssc} -rss {rss} -ng {ng}  -kfm {kfm} -nss {nss} -nsx {nsx} -nsy {nsy}'# -omi {omi}
	elif ngflag == True:
		#ng = f'3 6 {0.5*sigma_match} 4 {sigma_match} 2 {2.*sigma_match}'
		ng = '3 6 0.70 4 1.50 2 3.00'
		
		com = f'hotpants -c {conv} -n {scale} -iu {iu} -iuk {iu} -tu {tu} -tuk {tu} -il {il} -tl {tl} -v 0 -inim {sciim} -tmplim {refim} -outim {outfile} -oci {convfile} -ks {ks} -ko {ko} -bgo {bgo} -ft {ft} -ssig {ssig} -ig {igain} -tg {tgain} -ir {irdnoise} -tr {trdnoise} -afssc {afssc} -rss {rss} -ng {ng} -kfm {kfm} -nss {nss} -nsx {nsx} -nsy {nsy}'

	print(com)
	os.system(com)
	print('All done, check it out!')	
'''

def hotpants(sciim, refim, ngflag=False, sigma_match='', conv='t', scale='i', iu=5000000, tu=5000000, il=-10000, tl=-10000, igain=1.33, tgain=1.26, irdnoise=12.6, trdnoise=361., ssf='', **kwargs):
    """
    hotpants 실행 함수 (TARS)
    
    Parameters:
    sciim (str): 과학 이미지 파일명
    refim (str): 기준 이미지 파일명
    ngflag (bool): NG 설정 여부
    sigma_match (str): 매칭 시그마 값
    conv (str): 컨볼루션 옵션
    scale (str): 스케일 옵션
    iu, tu, il, tl (int): 픽셀 값 범위 설정
    igain, tgain (float): Gain 값 설정
    irdnoise, trdnoise (float): Read Noise 값 설정
    ssf (str): 기타 옵션
    **kwargs: 추가적인 설정값을 키워드 인자로 입력 가능

	# (1) Stamp, substamp

	nsx   = 12 # 이미지의 X축과 Y축에서 스탬프의 개수를 설정합니다. 기본값은 각각 10입니다.스탬프는 이미지의 특정 구역으로 나누어진 작은 영역이며, 각 영역에서 PSF와 커널을 개별적으로 계산합니다. (10) 
	nsy   = 12 # (10)

	nss   = 3 # 각 스탬프 내에서 사용할 서브스탬프 수를 설정합니다. 적절한 수를 설정해 문제 영역을 최소화할 수 있습니다. (3) 
	rss   = 15 # 각 스탬프의 반지름을 조정하여, 더 작은 영역을 검사하고 문제가 있는 부분을 피할 수 있습니다. (15) 


	# (2) Sigma clipping
	
	ks    = 2 # 커널 피팅 과정에서 이상치를 제거하기 위해 사용됩니다. 기본적으로 2.0으로 설정되며, 이 값은 커널 내 픽셀 값이 평균으로부터 몇 시그마 이상 벗어날 경우 이를 이상치로 간주하고 제거하는 기준입니다. (2)

	ssig  = 3 # -ks는 커널 피팅에서, -ssig는 통계적 분석에서 Sigma 클리핑 임계값을 설정하는 옵션입니다. 이 값을 통해 이상치나 노이즈가 있는 스탬프를 걸러낼 수 있으며, 이를 적절히 설정하지 않으면 NaN값이 발생할 수 있습니다. 이미지의 노이즈 특성에 맞게 최적의 값을 찾아야 합니다. (3.0) 

	# (3) Kernel order & Composition
	
	#r     = 10 # 컨볼루션 커널의 반지름을 설정하는 옵션입니다.커널의 크기는 이미지의 PSF 크기에 맞춰 조정해야 합니다. 너무 작은 값을 설정하면 PSF의 전체 범위를 포착하지 못하고, 너무 큰 값을 설정하면 불필요한 영역이 포함되어 노이즈가 증가할 수 있습니다.(10)

	ko    = 2 # 커널이 이미지 내에서 공간적으로 어떻게 변화하는지를 설명하는 다항식의 차수입니다.PSF가 이미지 내에서 변화할 때, 이 옵션을 적절히 설정하여 복잡한 PSF 변화를 반영할 수 있습니다. 잘못된 차수 설정은 PSF를 부정확하게 모델링할 수 있으므로, 상황에 맞게 조정하는 것이 중요합니다.(2) 

	bgo   = 1 # 이미지 내에서 배경광의 변화를 설명하는 다항식의 차수입니다. 배경이 균일하지 않은 이미지에서 적절한 값을 설정하면 배경 차이를 효과적으로 보정할 수 있습니다. 배경 변화를 정확하게 반영하지 않으면 차이 이미지에 잔여 패턴이나 거짓 신호가 남을 수 있습니다. (1) 

	kfm   = 0.99 #커널의 절대 값 합 중 유효한 값을 가질 수 있는 픽셀의 비율을 설정하는 옵션이에요. 기본값은 0.990로, 전체 커널의 99%를 유효한 값으로 간주합니다.

	# ng 

	# (4) RMS Threshold and Fit Threshold

	ft    =  21 # 커널 피팅에서 유효한 중심값을 결정하기 위한 RMS 임계값을 설정합니다. 기본값은 20으로 설정되어 있습니다. 이 값은 스탬프의 중심값(centroid)이 유효한지 여부를 결정하며, 특정 스탬프에서 이 임계값을 넘으면 해당 스탬프가 이상치로 간주될 수 있습니다. (20.0)	

	# ft = 7 # DOAO 필드 안에 별이 X, MAO_FLI도 별이 X
	afssc = 1 # 스탬프의 자동 중심 탐색을 활성화하거나 비활성화합니다.PSF가 잘 정의된 영역을 찾고, 스탬프의 중심을 자동으로 선택하여 PSF가 왜곡된 영역을 피할 수 있도록 설정하는 것이 중요합니다. (1)
    """
    # 기본값 설정
    params = {
        'nsx': 12, 'nsy': 12,  # Stamp 개수
        'nss': 3, 'rss': 15,    # Substamp 설정
        'ks': 2, 'ssig': 3,     # Sigma clipping
        'ko': 2, 'bgo': 1,      # Kernel order & Background order
        'kfm': 0.99,            # Kernel fraction mask
        'ft': 21, 'afssc': 1,   # RMS Threshold & Auto-fit stamp center
    }
    
    params.update(kwargs)
		
    outfile = 'hd' + sciim
    convfile = 'hc' + sciim
    
    # NG 설정
    if not ngflag:
        ng = '3 6 0.70 4 1.50 2 3.00'
    else:
        ng = '3 6 0.70 4 1.50 2 3.00'
    
    # hotpants 명령어 생성
    com = (f'hotpants -c {conv} -n {scale} -iu {iu} -iuk {iu} -tu {tu} -tuk {tu} -il {il} -tl {tl} '
           f'-v 0 -inim {sciim} -tmplim {refim} -outim {outfile} -oci {convfile} '
           f'-ks {params["ks"]} -ko {params["ko"]} -bgo {params["bgo"]} -ft {params["ft"]} -ssig {params["ssig"]} '
           f'-ig {igain} -tg {tgain} -ir {irdnoise} -tr {trdnoise} -afssc {params["afssc"]} -rss {params["rss"]} '
           f'-ng {ng} -kfm {params["kfm"]} -nss {params["nss"]} -nsx {params["nsx"]} -nsy {params["nsy"]}')
    
    print(com)
    os.system(com)
    print('All done, check it out!')


def runhotpants(imlist_name='Cal*fits', refinst='PS1', scifwhmkey='FWHM', reffwhmkey='FWHM', refim='', same=True, part0='Remap', **kwargs):
    imlist = glob.glob(imlist_name); imlist.sort()

    #obsspec = ascii.read('/home/lim9/miniconda3/lib/python3.9/site-packages/lgpyred/data/obs_spec.txt')
    with pkg_resources.path(lgpyred.data, 'obs_spec.txt') as path:
        obsspec = ascii.read( str(path))    
    
    refspec = obsspec[np.where(obsspec['obs_ccd'] == refinst)[0]]

    # 설정값 결정
    param_defaults = {
        'PS1': {'tu': 3000000, 'tl': -100000, 'il': -100, 'iu': 50000},
        'ZTF': {'tu': 50000, 'tl': 1, 'il': -100, 'iu': 50000},
        'LCOGT_BANZAI': {'tu': 50000, 'tl': -10, 'il': -100, 'iu': 50000},
        'LCOGT': {'tu': 50000, 'tl': -10, 'il': -100, 'iu': 50000},
        'LOAO': {'tu': 50000, 'tl': 1, 'il': -1000, 'iu': 30000},
        'LOAO_E2V': {'tu': 50000, 'tl': 1, 'il': -1000, 'iu': 30000},
        'DOAO': {'tu': 50000, 'tl': -100, 'il': 0, 'iu': 60000},
        'DOAO_SOPHIA': {'tu': 50000, 'tl': -100, 'il': 0, 'iu': 60000},
        'MAO_FLI': {'tu': 50000, 'tl': -500, 'il': -500, 'iu': 50000},
		'PNUO_C361K' : {'tu' : 50000, 'tl' : -100, 'il' : -100, 'iu' : 50000}	
    }

    # refinst가 param_defaults에 없으면 기본값을 사용
    params = param_defaults.get(refinst, {}).copy()

    # params가 None이 아니면 kwargs로 업데이트
    if params is not None:
        params.update(kwargs)
    else:
        params = kwargs  # refinst가 없으면 kwargs만 사용

    for inim in imlist:
        sciim = inim
        print(f'Science = {sciim}')
        
        part = sciim.split('-')
        sciinst = part[1]
        scispec = obsspec[np.where(obsspec['obs_ccd'] == sciinst)[0]]
        scihdr = fits.getheader(sciim, ignore_missing_end=True)
        scifwhm = scihdr[scifwhmkey]

        if same:
            pixscale = scispec['pixelscale'].item()
            part[0] = part0
            refim = '-'.join(part)
            print(f'Template = {refim}')
            refhdr = fits.getheader(refim, ignore_missing_end=True)
            reffwhm = refhdr[reffwhmkey]

            print(f'FWHM_ref = {reffwhm}"')
            print(f'FWHM_sci = {scifwhm}"')

            sigma_template = float(reffwhm) / 2.355
            sigma_image = float(scifwhm) / 2.355

            if sigma_image > sigma_template:
                sigma_match = np.sqrt((sigma_image / pixscale) ** 2 - (sigma_template / pixscale) ** 2)
            else:
                sigma_match = np.sqrt((sigma_template / pixscale) ** 2 - (sigma_image / pixscale) ** 2)

            stamp = sciim[:-8] + 'stamp.dat'

            if sigma_image >= sigma_template:
                print('Template image will be convolved.')
                hotpants(
                    sciim, refim, conv='t', scale='i', ngflag=True, sigma_match=sigma_match, igain=scispec['gain'].item(), tgain=refspec['gain'].item(), irdnoise=scispec['RDnoise'].item(), trdnoise=scispec['RDnoise'].item(), ssf=stamp,
                    **params  # 여기에 추가적인 설정 가능
                )

            elif float(scifwhm) < float(reffwhm):
                print('Science image will be convolved.')
                hotpants(
                    sciim, refim, conv='i', scale='i', ngflag=False, sigma_match='', igain=scispec['gain'].item(), tgain=refspec['gain'].item(),
                    irdnoise=scispec['RDnoise'].item(), trdnoise=scispec['RDnoise'].item(), ssf=stamp,
                    **params  # 여기에 추가적인 설정 가능
                )

        elif not same:
            pixscale = refspec['pixelscale'].item()
            print(f'Template = {refim}')
            refhdr = fits.getheader(refim, ignore_missing_end=True)
            reffwhm = refhdr[reffwhmkey]

            print(f'FWHM_ref = {reffwhm}"')
            print(f'FWHM_sci = {scifwhm}"')

            sigma_template = float(reffwhm) / 2.355
            sigma_image = float(scifwhm) / 2.355

            if sigma_image > sigma_template:
                sigma_match = np.sqrt((sigma_image / pixscale) ** 2 - (sigma_template / pixscale) ** 2)
            else:
                sigma_match = np.sqrt((sigma_template / pixscale) ** 2 - (sigma_image / pixscale) ** 2)

            stamp = inim[:-8] + 'stamp.dat'

            if sigma_image >= sigma_template:
                print('Template image will be convolved.')
                hotpants(
                    sciim, refim, conv='t', scale='i', ngflag=True, sigma_match=sigma_match, igain=scispec['gain'].item(), tgain=refspec['gain'].item(),
                    irdnoise=scispec['RDnoise'].item(), trdnoise=scispec['RDnoise'].item(), ssf=stamp,
                    **params  # 여기에 추가적인 설정 가능
                )

            elif sigma_image < sigma_template:
                print('Science image will be convolved.')
                hotpants(
                    sciim, refim, conv='i', scale='i', ngflag=False, sigma_match='', igain=scispec['gain'].item(), tgain=refspec['gain'].item(),
                    irdnoise=scispec['RDnoise'].item(), trdnoise=scispec['RDnoise'].item(), ssf=stamp,
                    **params  # 여기에 추가적인 설정 가능
                )

    print('All images are subtracted.')

'''
def runhotpants(imlist_name='Cal*fits', refinst = 'PS1', scifwhmkey='FWHM', reffwhmkey='FWHM', refim='', same=True, part0='Remap') :

	imlist = glob.glob(imlist_name); imlist.sort()

	obsspec = ascii.read('/home/lim9/miniconda3/lib/python3.9/site-packages/lgpyred/data/obs_spec.txt')

	refspec = obsspec[ np.where( obsspec['obs_ccd'] ==  refinst )[0]]

	for inim in imlist:
		sciim  = inim
		print(f'Science = {sciim}')
		
		part    = sciim.split('-')

		sciinst = part[1]
		scispec = obsspec[ np.where( obsspec['obs_ccd'] ==  sciinst )[0]]
		scihdr  = fits.getheader(sciim, ignore_missing_end=True)
		scifwhm = scihdr[scifwhmkey]


		if refinst == 'PS1':
			tu = 3000000
			tl = -100000
			il = -100
			iu = 50000

		elif refinst == 'ZTF' :
			tu = 50000
			tl = 1
			il = -100
			iu = 50000
		elif (refinst == 'LCOGT_BANZAI') | (refinst == 'LCOGT') :
			tu = 50000.
			tl = -10		
			il = -100
			iu = 50000	

		elif (refinst == 'LOAO') | (refinst == 'LOAO_E2V') :
			tu = 50000
			tl = 1	
			il = -1000 # LOAO I with fringe pattern
			iu = 30000
		elif (refinst == 'DOAO') | (refinst == 'DOAO_SOPHIA') :
			tu = 50000
			tl = -100	
			il = 0
			iu = 60000			
		elif refinst == 'MAO_FLI' :
			tu = 50000
			tl = -500	
			il = -500
			iu = 50000		
		else:
			tu = 50000
			tl = -100	
			il = -100
			iu = 50000		

		if same == True:
			pixscale = scispec['pixelscale'].item()
			
			part[0] = part0

			refim   = '-'.join(part)
			print(f'Template = {refim}')
			refhdr  = fits.getheader(refim, ignore_missing_end=True)
			reffwhm = refhdr[reffwhmkey]

			print(f'FWHM_ref = {reffwhm}"')
			print(f'FWHM_sci = {scifwhm}"')

			stamp = sciim[:-8] + 'stamp.dat'

			sigma_template = float(reffwhm)/2.355
			print(f'sigma (template) = {sigma_template:.2f}"')
			sigma_image    = float(scifwhm)/2.355
			print(f'sigma (image) = {sigma_image:.2f}"')
			if sigma_image > sigma_template :
				sigma_match = np.sqrt((sigma_image/pixscale)**2 - (sigma_template/pixscale)**2)
			elif sigma_image <= sigma_template :
				sigma_match = np.sqrt((sigma_template/pixscale)**2 - (sigma_image/pixscale)**2)

			if sigma_image >= sigma_template :
				print('Template image will be convolved.')
				# This leads to smoothing of the template. Assume that both Psfs are Gaussian, in which case the Gaussian that matches the two has Sigma_match = sqrt(Sigma_image2 - Sigma_template2). It is recommended that this be the central Gaussian in your kernel basis, with the smallest one being 0.5 * Sigma_match and the largest being 2.0 * Sigma_match. Set these using the -ng flag. E.g. -ng 3 6 0.5Sigma_match 4 Sigma_match 2 2.0Sigma_match.

				hotpants(sciim, refim, conv='t', scale='i', ngflag=True, sigma_match=sigma_match, iu=iu, tu=tu, il=il, tl=tl, igain=scispec['gain'].item(), tgain=refspec['gain'].item(), irdnoise=scispec['RDnoise'].item(), trdnoise=scispec['RDnoise'].item() , ssf=stamp) 

			elif float(scifwhm) < float(reffwhm) :
				print('Science image will be convolved.')
				hotpants(sciim, refim, conv='i', scale='i', ngflag=False, sigma_match='', iu=iu, tu=tu, il=il, tl=tl, igain=scispec['gain'].item(), tgain=refspec['gain'].item(), irdnoise=scispec['RDnoise'].item(), trdnoise=scispec['RDnoise'].item(),  ssf=stamp)
		elif same == False:
			pixscale = refspec['pixelscale'].item()
			print(f'Template = {refim}')
			refhdr  = fits.getheader(refim, ignore_missing_end=True)
			reffwhm = refhdr[reffwhmkey]

			print(f'FWHM_ref = {reffwhm}"')
			print(f'FWHM_sci = {scifwhm}"')
			sigma_template = float(reffwhm)/2.355
			print(f'sigma (template) = {sigma_template:.2f}"')
			sigma_image    = float(scifwhm)/2.355
			print(f'sigma (image) = {sigma_image:.2f}"')

			if sigma_image > sigma_template :
				sigma_match = np.sqrt((sigma_image/pixscale)**2 - (sigma_template/pixscale)**2)
			elif sigma_image <= sigma_template :
				sigma_match = np.sqrt((sigma_template/pixscale)**2 - (sigma_image/pixscale)**2)

			stamp = inim[:-8] + 'stamp.dat'

			if sigma_image >= sigma_template :
				print('Template image will be convolved.')
				hotpants(sciim, refim, conv='t', scale='i', ngflag=True, sigma_match=sigma_match, iu=iu, tu=tu, il=il, tl=tl, igain=scispec['gain'].item(), tgain=refspec['gain'].item(), irdnoise=scispec['RDnoise'].item(), trdnoise=scispec['RDnoise'].item(),  ssf=stamp) 

			elif sigma_image < sigma_template :
				print('Science image will be convolved.')
				hotpants(sciim, refim, conv='i', scale='i', ngflag=False, sigma_match='', iu=iu, tu=tu, il=il, tl=tl, igain=scispec['gain'].item(), tgain=refspec['gain'].item(), irdnoise=scispec['RDnoise'].item(), trdnoise=scispec['RDnoise'].item(),  ssf=stamp)
	print('All images are subtracted.')
'''
