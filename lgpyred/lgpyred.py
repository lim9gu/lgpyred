### lgpyobs image reduction pipeline

import time
import shutil
import os, sys, glob
import argparse
import datetime
from itertools import product
from datetime import datetime as dt

import paramiko
import numpy as np
from pyraf import iraf
import astropy.units as u
from astropy.time import Time
from astropy.table import Table
from astropy.table import Column
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord

from lgpytars.SameEpoch import SameEpoch
from lgpytars.imcombine import imcombine_set, imcombine_epoch
from lgpytars.reduction.hdrcheck import wcscenter
from lgpytars.reduction.wcsremap import run_wcsremap_epoch, run_wcsremap_wrap
from lgpytars.lgphot import phot
from lgpytars.reduction import hotpants

class Red :
    def __init__(self,
                 imlist_name = '*.fit',
                 sumfile     = 'file_summary.txt',
                 ccd         = 'PNUO_C361K',
                 zeroproc    = True, 
                 darkproc    = True,
                 flatproc    = True,
                 flattype    = 'skyflat',
                 sep         = 60.
                 ) :
        """
        Initialize the Red class with settings and configurations for image reduction.

        Parameters:
        imlist_name (str): Pattern to match image files.
        zeroproc (bool): Flag to indicate if zero (bias) processing should be done.
        darkproc (bool): Flag to indicate if dark processing should be done.
        flattype (str): Type of flat field ('skyflat' or 'Domeflat').

        Raises:
        ValueError: If the flattype is not 'skyflat' or 'Domeflat'.
        """
        # File setting
        self.curdir      = os.getcwd()

        iraf.chdir(self.curdir)
        self.irafdb = '/home/lim9/iraf-2.17/noao/imred/ccdred/ccddb/'
        self.imlist_name = imlist_name
        if '-' in list(self.curdir.split('/')[-1]) :
            yy   = self.curdir.split('/')[-1].split('-')[0]
            mm   = self.curdir.split('/')[-1].split('-')[1]
            dd   = self.curdir.split('/')[-1].split('-')[2]
            self.curdate = yy+mm+dd
        else :
            self.curdate = self.curdir.split('/')[-1]

        # Detector setting
        self.ccd         = ccd
        obscat           = ascii.read("/home/lim9/miniconda3/lib/python3.9/site-packages/lgpytars/data/obs_spec.txt")
        inst             = obscat[obscat['obs_ccd'] == self.ccd]

        self.temp2replace = -5. # If no CCD-TEMP, replace this temp.

        self.gain        = inst['gain'][0]
        self.rdnoise     = inst['RDnoise'][0]
        self.pixscale    = inst['pixelscale'][0]
        self.magl        = inst['magl'][0]
        self.magh        = inst['magh'][0]
        self.fov         = inst['fov_a'][0]
        self.savedir     = inst['savedir'][0]
        self.sumfile     = sumfile
        #self.calfile     = 'calib_summary.txt'

        # Reduction setting
        self.zeroproc    = True
        self.darkproc    = True
        self.flatproc    = True

        self.flat_zeroproc = True
        self.flat_darkproc = True

        if flattype not in ['skyflat', 'Domeflat'] :
            raise ValueError(f"Invalid flattype: {flattype}. Must be 'skyflat' or 'Domeflat'.")
        self.flattype    = flattype

        self.bias_archive   = f'{self.savedir}masterbias/'
        self.dark_archive   = f'{self.savedir}masterdark/'
        self.flat_archive   = f'{self.savedir}masterflat_*/'

        # Astrometry setting (Astrometry.net)
        self.config    = '/home/lim9/miniconda3/etc/astrometry.cfg'
        self.scalelow  = self.pixscale * 0.9 #0.18  # lower limit for the pixel scale covering 'No binning' to '2x2 binning')
        self.scalehigh = self.pixscale * 2 * 1.1 #0.46  # upper limit for the pixel scale covering 'No binning' to '2x2 binning')
        self.radec     = False # if True, astronometry.net search for wcs using radec in fits image header 
        self.sexconfig = f'/home/lim9/miniconda3/lib/python3.9/site-packages/lgpytars/astrom.config/astrometry.net.{self.ccd}.sex'
        self.radius    = 0.7 # degree of querying radius of stars in index files

        # Header check
        self.alltarget = "/home/lim9/miniconda3/lib/python3.9/site-packages/lgpytars/data/alltarget.dat"

        # Image alignment setting (wcsremap)
        self.sep         = sep   # mins time interval
        self.wrapper     = True # use wcsremap wrapper?

        # Image stacking setting (IRAF.imcombine)
        self.reject = 'none' 

        # Image subtraction setting (HOTPANTS)
        self.template_dir = '/mnt/dataset/obsdata/IMSNG/template_20250213'

        # File transfer setting
        self.server = ''

        self.banner()

        print(f'Reduction for {self.curdir} is now ready.')

    #@staticmethod
    def banner(self):
        print(rf"""
        ██╗      ██████╗ ██████╗ ██╗   ██╗██████╗ ███████╗██████╗ 
        ██║     ██╔════╝ ██╔══██╗╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
        ██║     ██║  ███╗██████╔╝ ╚████╔╝ ██████╔╝█████╗  ██║  ██║
        ██║     ██║   ██║██╔═══╝   ╚██╔╝  ██╔══██╗██╔══╝  ██║  ██║
        ███████╗╚██████╔╝██║        ██║   ██║  ██║███████╗██████╔╝
        ╚══════╝ ╚═════╝ ╚═╝        ╚═╝   ╚═╝  ╚═╝╚══════╝╚═════╝ 
        """)
        print(rf"""
        LGPY imaging REDuction pipeline for {self.ccd} 
        2025-01-02 Version (by G. Lim; lim9gu@gmail.com) 
        """)
        # banner from https://patorjk.com/software/taag/#p=display&v=0&f=ANSI%20Shadow&t=lgpyred%0A-PNUO

    def GetMaster(self, filetype, bin=None, exp=None, temp=None, filter_name=None):
        """
        Find the closest master calibration file (bias, dark, or flat) based on the date.

        Parameters:
        - filetype: Type of the master file ('bias', 'dark', 'flat')
        - bin: Binning setting (optional for dark/flat)
        - exp: Exposure time (optional for dark/flat)
        - temp: CCD temperature (optional for dark)

        Returns:
        - The path of the closest master calibration file
        """


        # 파일 경로 패턴 정의
        if filetype == 'flat' and filter_name is not None:
            pattern = f"{self.savedir}masterflat_{filter_name}/*_n{filter_name}skyflat_bin{bin}.fits"
        elif (filetype == 'bias') | (filetype == 'zero'):
            pattern = f"{self.savedir}masterbias/*_zero_bin{bin}.fits"
        elif filetype == 'dark':
            pattern = f"{self.savedir}masterdark/*_dark{exp}_bin{bin}_{temp}.fits"
        else:
            raise ValueError("Invalid filetype or missing parameters for master flat.")

        # 파일 검색
        master_files = glob.glob(pattern)
        if not master_files:
            print(f"No master {filetype} files found for filter: {filter_name}.")
            return None

        # 날짜 차이 계산
        date_diffs = []
        today = dt.strptime(self.curdate, '%Y%m%d')

        for file in master_files:
            try:
                # 파일 이름에서 날짜 추출
                date_str = os.path.basename(file).split('_')[0]
                file_date = dt.strptime(date_str, '%Y%m%d')
                diff = abs((file_date - today).days)  # 날짜 차이 계산
                date_diffs.append((diff, file))
            except Exception as e:
                print(f"Error parsing date from file: {file}, {e}")
                continue

        # 가장 가까운 날짜 파일 선택
        if date_diffs:
            date_diffs.sort(key=lambda x: x[0])  # 차이로 정렬
            closest_file = date_diffs[0][1]     # 가장 가까운 파일 반환
            print(f"Closest master {filetype} found: {closest_file}")
            return closest_file
        else:
            print(f"No valid master {filetype} files with dates found.")
            return None

    def FileSum(self, imlist_name=None) :
        """
        Create a summary of the image files in the current directory.

        This method reads image files matching `imlist_name`, extracts relevant
        header information, and creates a summary table which is saved to `self.sumfile`.

        The summary includes:
        - Filename
        - Modified Julian Date (MJD)
        - X and Y binning
        - Image type (Bias, Dark, Flat, Light)
        - Exposure time
        - Filter used
        - Object name

        Additionally, it updates the image headers to replace spaces in the object name
        and adds the MJD.
        """
        #imlist = glob.glob(self.imlist_name); imlist.sort()
        if self.imlist_name == None:
            imlist = self.imlist
        elif self.imlist_name != None:
            imlist = glob.glob(self.imlist_name); imlist.sort()

        XBIN, YBIN, FILENAME, IMAGETYP, EXPTIME = [], [], [], [], []
        FILTER, OBJECT, MJD, CCDTEMP = [], [], [], []

        # Making image info table
        for im in imlist :
            print(im)

            # Remove space in name
            newim = im.replace(" ","")
            if im != newim :
                print(f'{im} -> {newim}')
                os.rename(im, newim)

            hdr   = fits.getheader(newim)

            XBIN.append(hdr['XBINNING'])
            YBIN.append(hdr['YBINNING'])
            FILENAME.append(newim)
            IMAGETYP.append(hdr['IMAGETYP'])
            EXPTIME.append(hdr['EXPTIME'])

            if hdr['IMAGETYP'] == 'Bias Frame' :
                FILTER.append('Bias')
            elif hdr['IMAGETYP'] == 'Dark Frame' :
                FILTER.append('Dark')
            elif hdr['IMAGETYP'] == 'Flat Field' :
                FILTER.append(hdr['FILTER'])
            elif hdr['IMAGETYP'] == 'Light Frame' :
                FILTER.append(hdr['FILTER'])

            if hdr['OBJECT'] == "" or hdr['OBJECT'] is None:
                hdr['OBJECT'] = 'NONAME'
            
            if 'CCD-TEMP' not in hdr :
                print(f'{im} does not have CCD-TEMP in header.')
                hdr['CCD-TEMP'] = self.temp2replace

            OBJECT.append(hdr['OBJECT'].replace(" ","").replace("-",""))

            CCDTEMP.append(round(hdr['CCD-TEMP'], 2))

            t     = Time(hdr['DATE-OBS'], format='isot', scale='utc')
            MJD.append(round(t.mjd,6))
            
            iraf.hedit(newim, fields='OBJECT', value=hdr['OBJECT'].replace(" ","").replace("-",""), verify='no')
            iraf.hedit(newim, fields='MJD', value=round(t.mjd,6), verify='no', add='yes')
            
        cat = Table({
            'FILENAME': FILENAME,
            'MJD': MJD,
            'XBINNING': XBIN,
            'YBINNING': YBIN,
            'IMAGETYP': IMAGETYP,
            'EXPTIME': EXPTIME,
            'FILTER': FILTER,
            'OBJECT': OBJECT,
            'CCDTEMP': CCDTEMP,
        }, names=['FILENAME', 'MJD', 'XBINNING', 'YBINNING', 'IMAGETYP', 'EXPTIME', 'FILTER', 'OBJECT', 'CCDTEMP'])

        cat['CCDTEMP0'] = np.round( cat['CCDTEMP'] ).astype(int)

        ascii.write(cat, self.sumfile, format='basic', overwrite=True)
        print("File summary is created.")

        # Finding missing binning, exptime, filter.
        bias = cat[cat['IMAGETYP'] == 'Bias Frame']
        dark = cat[cat['IMAGETYP'] == 'Dark Frame']
        flat = cat[cat['IMAGETYP'] == 'Flat Field']
        sci  = cat[cat['IMAGETYP'] == 'Light Frame']

        for obj in sci:
            binning, exptime, band, ccdtemp = obj['XBINNING'], obj['EXPTIME'], obj['FILTER'], obj['CCDTEMP0']

            # Check for raw data first
            raw_bias = bias[(bias['XBINNING'] == binning)]
            raw_dark = dark[(dark['XBINNING'] == binning) & (dark['EXPTIME'] == exptime) & (dark['CCDTEMP0'] == ccdtemp)]
            raw_flat = flat[(flat['XBINNING'] == binning) & (flat['FILTER'] == band)]

            # Handle Bias
            if len(raw_bias) > 0:
                print(f"Raw bias exists for bin: {binning}.")
            else:
                print(f"No raw bias found for bin: {binning}.")
                closest_bias = self.GetMaster('bias', bin=binning)
                if closest_bias:
                    print(f"Using closest bias: {closest_bias}")
                    os.system(f'cp {closest_bias} ./')
                else:
                    print("No bias available.")

            # Handle Dark
            if len(raw_dark) > 0:
                print(f"Raw dark exists for bin: {binning}, exp: {exptime}, temp: {ccdtemp}.")
            else:
                print(f"No raw dark found for bin: {binning}, exp: {exptime}, temp: {ccdtemp}.")
                closest_dark = self.GetMaster('dark', bin=binning, exp=exptime, temp=ccdtemp)
                if closest_dark:
                    print(f"Using closest dark: {closest_dark}")
                    os.system(f'cp {closest_dark} ./')
                else:
                    print("No dark available.")

            # Handle Flat
            if len(raw_flat) > 0:
                print(f"Raw flat exists for bin: {binning}, band: {band}.")
            else:
                print(f"No raw flat found for bin: {binning}, band: {band}.")
                closest_flat = self.GetMaster('flat', bin=binning, filter_name=band)
                if closest_flat:
                    print(f"Using closest flat: {closest_flat}")
                    os.system(f'cp {closest_flat} ./')
                else:
                    print("No flat available.")

        '''
        print('Inspecting dark...')
        for obj in dark :
            binning, exptime, ccdtemp = obj['XBINNING'], obj['EXPTIME'], obj['CCDTEMP0']
            #exptime = obj['EXPTIME']
            #band    = obj['FILTER']
            closest_bias = self.GetMaster('bias', bin=binning)
            if closest_bias:
                print(f"Closest bias file found: {closest_bias}")
                os.system(f'cp {closest_bias} ./')
            else:
                print(f"No bias found for bin: {binning}.") 
            #"""       
            if not ((bias['XBINNING'] == binning)).any():
                #missing_bin.append(obj['XBINNING'])

                print(f'No bias (bin: {binning}) exists.')
                mbias = glob.glob(f'{self.bias_archive}2*_zero_bin{binning}.fits'); mbias.sort()
                try:
                    os.system(f'cp {mbias[-1]} ./')
                except:
                    print('No bias in archive.')
            else:
                pass
            #"""
        print('Inspecting flat...')
        for obj in flat :
            binning, band, ccdtemp = obj['XBINNING'], obj['FILTER'], obj['CCDTEMP0']
            closest_bias = self.GetMaster('bias', bin=binning)
            closest_dark = self.GetMaster('dark', bin=binning, exp=obj['EXPTIME'], temp=ccdtemp)
            closest_flat = self.GetMaster('flat', bin=binning, filter_name=band)

            if closest_bias:
                print(f"Closest bias file found: {closest_bias}")
                os.system(f'cp {closest_bias} ./')
            else:
                print(f"No bias found for bin: {binning}.")

            if closest_dark:
                print(f"Closest dark file found: {closest_dark}")
                os.system(f'cp {closest_dark} ./')
            else:
                print(f"No dark found for bin: {binning}, exp: {obj['EXPTIME']}, temp: {ccdtemp}.")

            if closest_flat:
                print(f"Closest flat file found: {closest_flat}")
                os.system(f'cp {closest_flat} ./')
            else:
                print(f"No flat found for bin: {binning}, band: {band}.")

            #"""
            binning = obj['XBINNING']
            exptime = obj['EXPTIME']
            band    = obj['FILTER']
            ccdtemp = obj['CCDTEMP0']

            if not ((bias['XBINNING'] == binning)).any() :
                #missing_bin.append(obj['XBINNING'])
                print(f'No bias (bin: {binning}) exists.')
                mbias = glob.glob(f'{self.bias_archive}2*_zero_bin{binning}.fits'); mbias.sort()
                try:
                    os.system(f'cp {mbias[-1]} ./')
                except:
                    print('No bias in archive.')
            else:
                pass
                
            if not ((dark['XBINNING'] == binning) & (dark['EXPTIME'] == exptime) & (dark['CCDTEMP0'] == ccdtemp)).any() :
                print(f'No dark (bin: {binning}, exp: {exptime}, ccdtemp: {ccdtemp}) exists.')
                mdark = glob.glob(f'{self.dark_archive}2*_dark{exptime}_bin{binning}_{ccdtemp}.fits'); mdark.sort()
                try:
                    os.system(f'cp {mdark[-1]} ./')
                except:
                    print('No dark in archive.')
            else:
                pass
            #"""

        print('Inspecting science images...')
        for obj in sci :
            binning, exptime, band, ccdtemp = obj['XBINNING'], obj['EXPTIME'], obj['FILTER'], obj['CCDTEMP0']
            closest_bias = self.GetMaster('bias', bin=binning)
            closest_dark = self.GetMaster('dark', bin=binning, exp=exptime, temp=ccdtemp)
            closest_flat = self.GetMaster('flat', bin=binning, filter_name=band)

            if closest_bias:
                print(f"Closest bias file found: {closest_bias}")
                os.system(f'cp {closest_bias} ./')
            else:
                print(f"No bias found for bin: {binning}.")

            if closest_dark:
                print(f"Closest dark file found: {closest_dark}")
                os.system(f'cp {closest_dark} ./')
            else:
                print(f"No dark found for bin: {binning}, exp: {exptime}, temp: {ccdtemp}.")

            if closest_flat:
                print(f"Closest flat file found: {closest_flat}")
                os.system(f'cp {closest_flat} ./')
            else:
                print(f"No flat found for bin: {binning}, band: {band}.")

            #"""
            binning = obj['XBINNING']
            exptime = obj['EXPTIME']
            band    = obj['FILTER']
            ccdtemp = obj['CCDTEMP0']

            if not ((bias['XBINNING'] == binning) ).any():
                #missing_bin.append(obj['XBINNING'])
                print(f'No bias (bin: {binning}) exists.')
                mbias = glob.glob(f'{self.bias_archive}2*_zero_bin{binning}.fits'); mbias.sort()
                try:
                    os.system(f'/usr/bin/cp {mbias[-1]} ./')
                except:
                    print('No bias in archive.')
            else:
                pass

            if not ((dark['XBINNING'] == binning) & (dark['EXPTIME'] == exptime) & (dark['CCDTEMP0'] == ccdtemp)).any() :
                #missing_bin.append(obj['XBINNING'])
                #missing_exp.append(obj['EXPTIME'])
                print(f'No dark (bin: {binning}, exp: {exptime}, ccdtemp: {ccdtemp}) exists.')
                mdark = glob.glob(f'{self.dark_archive}2*_dark{exptime}_bin{binning}_{ccdtemp}.fits'); mdark.sort()
                try:
                    os.system(f'/usr/bin/cp {mdark[-1]} ./')
                except:
                    print('No dark in archive.')
            else:
                pass

            if not ((flat['XBINNING'] == binning) & (flat['FILTER'] == band) ).any() :
                #missing_bin.append(obj['XBINNING'])
                #missing_exp.append(obj['EXPTIME'])
                #missing_band.append(obj['FILTER'])
                print(f'No flat (bin: {binning}, exp: {exptime}, band: {band}) exists.')
                mflat = glob.glob(f'{self.flat_archive}2*_n{band}skyflat_bin{binning}.fits'); mflat.sort()
                try:
                    os.system(f'/usr/bin/cp {mflat[-1]} ./')
                except:
                    print('No flat in archive.')
            else:
                pass
            #"""
        '''
        print('Inspection finished.')

        self.summary = cat
        self.summary.sort('MJD')
    
    def GenBias(self):
        """
        Generate a master bias frame from individual bias frames.

        This method:
        - Reads the summary table to identify bias frames.
        - Groups bias frames by their binning settings.
        - Combines the frames using median combination to create a master bias frame.
        - Saves the master bias frame to the specified directory.

        It uses IRAF tasks to perform image statistics and combination.

        Notes:
        - The master bias frame is saved as `masterbias/{curdate}_zero_bin{bin}.fits`.
        - This method removes temporary list files after processing.
        """
        iraf.noao(); iraf.imred(); iraf.ccdred()
        iraf.ccdred.setinst(instrume='camera', directo=self.irafdb, query='q', review='no')

        sum     = self.summary
        curdate = self.curdate
        savedir = self.savedir

        bias    = sum[sum['IMAGETYP'] == 'Bias Frame']
        binset  = list(set(bias['XBINNING'])) ; binset.sort()

        for bin in binset : 
            biasset   = list(set( bias['FILENAME'][ (bias['IMAGETYP'] == 'Bias Frame' ) & (bias['XBINNING'] == bin)] )) ; biasset.sort()

            if biasset == [] :
                continue

            input_name  = f'bias{bin}.list'
            output_name = f'{curdate}_zero_bin{bin}.fits'

            f = open(input_name,'w+')
            for b in biasset:
                f.write(b+'\n')
            f.close()      

            print('Zerocombine is running...')
            iraf.imstat(images='@'+input_name, fields='image,mean,midpt,mode,stddev,min,max', lsigma=3, usigma=3)
            if not os.path.exists(f'{savedir}{output_name}'):
                iraf.zerocombine(input='@'+input_name, output=output_name, combine='median', reject='none', process='no', scale='none', ccdtype='')
            else:
                print(f'{savedir}{output_name} already exists. Remove.')
                os.remove(f'{curdir}/{output_name}')
                iraf.zerocombine(input='@'+input_name, output=output_name, combine='median', reject='none', process='no', scale='none', ccdtype='')
            print(f'Output master {output_name} is created.')
            if not os.path.exists(f'{savedir}masterbias') :
                os.mkdir(f'{savedir}masterbias')
            shutil.copyfile(f'{self.curdir}/{output_name}', f'{savedir}masterbias/{output_name}')
        os.system('/usr/bin/rm bias*.list')
        iraf.dir('.')
        #self.masterbias = output_name

    def GenDark(self) :
        """
        Generate master dark frames from individual dark frames.

        This method:
        - Reads the summary table to identify dark frames.
        - Groups dark frames by their binning settings and exposure times.
        - Combines frames with the same binning and exposure time using median combination
          to create master dark frames.
        - Saves the master dark frames to the specified directory.

        Notes:
        - The master dark frames are saved with names indicating their exposure times, e.g.,
          `masterdark/{curdate}_dark{exp}.fits`.
        - This method removes temporary list files after processing.
        """

        sum     = self.summary
        curdate = self.curdate
        savedir = self.savedir

        dark       = sum[sum['IMAGETYP'] == 'Dark Frame']
        binset     = list(set(dark['XBINNING'])) ; binset.sort()
        darkexpset = list(set( dark['EXPTIME'] )); darkexpset.sort()
        darktempset =  list(set( dark['CCDTEMP0'] )); darktempset.sort()

        masterdark = []

        for bin, exp, temp in product(binset, darkexpset, darktempset):

            bias_file = self.GetMaster('bias', bin=bin)
            if not bias_file:
                print(f"No valid bias file found for bin: {bin}. Skipping...")
                continue
            darkset = sorted(set(dark['FILENAME'][(dark['IMAGETYP'] == 'Dark Frame') & (dark['XBINNING'] == bin) & (dark['EXPTIME'] == exp) &  (dark['CCDTEMP0'] == temp)]))

            if not darkset:
                print(f"No raw dark frames found for bin: {bin}, exp: {exp}, temp: {temp}. Skipping...")
                continue

            # Zero subtraction
            print('zero subtraction...')
            zdarkset = []
            for d in darkset:
                #iraf.imarith(operand1=d, op='-', operand2=f"{self.curdir}/{mblist[-1].split('/')[-1]}", result='z' + d)
                zoutput = 'z' + d
                iraf.imarith(operand1=d, op='-', operand2=bias_file, result=zoutput)
                zdarkset.append('z' + d)

            input_name = f'dark{exp}_{temp}.list'
            output_name = f'{curdate}_dark{exp}_bin{bin}_{temp}.fits'
            masterdark.append(output_name)

            with open(input_name, 'w+') as f:
                for zd in zdarkset:
                    f.write(zd + '\n')

            # Dark combine
            print('Darkcombine is running...')
            iraf.imstat(images='@'+input_name, fields='image,mean,midpt,mode,stddev,min,max', lsigma=3, usigma=3)

            if not os.path.exists(f'{savedir}{output_name}'):
                iraf.darkcombine(input='@'+input_name, output=output_name, combine='median', reject='crreject', process='no', scale='none', ccdtype='')
            else:
                print(f'{savedir}{output_name} already exists. Remove.')
                os.remove(f'{curdir}/{output_name}')
                iraf.darkcombine(input='@'+input_name, output=output_name, combine='median', reject='crreject', process='no', scale='none', ccdtype='')

            print(f'Output master {output_name} is created.')
            if not os.path.exists(f'{savedir}masterdark') :
                os.mkdir(f'{savedir}masterdark')
            shutil.copyfile(f'{self.curdir}/{output_name}', f'{savedir}masterdark/{output_name}')

        os.system('/usr/bin/rm d*.list')  
        os.system('/usr/bin/rm zcal*d*.fit')    
        iraf.imstat(images='2*d*.fits', fields='image,mean,midpt,mode,stddev,min,max', lsigma=3, usigma=3, Stdout=0)  
        iraf.dir('.')
        self.masterdark = masterdark

    def GenFlat(self) :
        """
        Generate master flat field images by processing raw flat field images.

        This function performs bias and dark subtraction on raw flat field images,
        combines them to create master flat fields, and normalizes these master
        flats.

        Steps involved:
        1. Filter the summary table for 'Flat Field' images.
        2. Extract unique binning settings, filter bands, and flat field object types.
        3. Perform bias and dark subtraction on each flat field image.
        4. Combine processed flat field images into master flats for each combination of binning, filter, and object type.
        5. Normalize the combined master flats.
        6. Save the normalized master flats to a specified directory.

        Attributes:
        - `summary`: A DataFrame containing metadata of images.
        - `curdate`: Current date used in naming the output files.
        - `savedir`: Directory where the master flat files will be saved.
        - `masterflat`: List of normalized master flat files created.

        Parameters:
        None

        Returns:
        None

        Raises:
        IOError: If there are issues with reading or writing image files.
        """
        sum       = self.summary
        curdate   = self.curdate
        savedir   = self.savedir

        flat      = sum[sum['IMAGETYP'] == 'Flat Field']

        flatbinset    = list(set(flat['XBINNING'])) ; flatbinset.sort()
        flatbandset   = list(set(flat['FILTER'][flat['IMAGETYP'] == 'Flat Field' ] ))
        flatobjectset = list(set(flat['OBJECT'])) # skyflat or Domeflat

        ### Bias, Dark subtraction
        for im in flat['FILENAME'] :
            print(im)
            imidx  = np.where(flat['FILENAME'] == im)[0]
            bin    = flat['XBINNING'][imidx][0]
            exp    = flat['EXPTIME'][imidx][0]
            band   = flat['FILTER'][imidx][0]
            ccdtemp = flat['CCDTEMP0'][imidx][0]

            bias_file = self.GetMaster('bias', bin=bin)
            dark_file = self.GetMaster('dark', bin=bin, exp=exp, temp=ccdtemp)            
            if not bias_file:
                raise IOError(f"No valid bias file found for bin: {bin}.")
            if not dark_file:
                raise IOError(f"No valid dark file found for bin: {bin}, exp: {exp}, temp: {ccdtemp}.")

            if self.flat_zeroproc :
                print('Zero subtraction...')
                #iraf.imarith(operand1=im, op='-', operand2=f"{self.curdir}/{mblist[-1].split('/')[-1]}", result='z'+im)
                iraf.imarith(operand1=im, op='-', operand2=bias_file, result='z' + im)
            else :
                print('Zero subtraction is skipped...')
                os.system(f'cp {im} z{im}')

            if self.flat_darkproc :
                print('Dark subtraction...')
                zinput = 'z'+im
                iraf.imarith(operand1=zinput, op='-', operand2=dark_file, result='d' + zinput)

            else :
                print('Dark subtraction is skipped...')
                os.system(f'cp z{im} dz{im}')

        ### Flat combine
        for bin in flatbinset :
            for band in flatbandset :
                for self.flattype in flatobjectset :
                    dzflat = iraf.hselect('dz*.fit', fields="$I", expr=f'OBJECT = "{self.flattype}" && IMAGETYP = "Flat Field" && FILTER = "{band}" && XBINNING = {bin}', Stdout = 1)

                    if dzflat == [] :
                        pass
                    else :
                        input_name  = f'{self.flattype}{band}.list'

                        f = open(input_name,'w+')    
                        for im in dzflat : 
                            f.write(im+'\n')
                        f.close()  

                        # Flat combine
                        print(f'input : {input_name}')
                        output_name = f'{input_name[:-5]}_bin{bin}.fits'
                        #masterflat.append(output_name)

                        print(f'output : {output_name}')
                        iraf.imstat(images='@'+input_name,fields='image,mean,midpt,mode,stddev,min,max', lsigma=3, usigma=3)
                        iraf.flatcombine(input='@'+input_name, output=output_name, combine='median', reject='avsigclip', process='no', scale='mode', ccdtype='', lsigma='3.', hsigma='3.')

                        print(f'{output_name} is created. Normalizing...')
                        data, newhdr = fits.getdata(output_name, header=True)
                        #x           = np.mean(data)
                        x            = np.median(data)
                        nimage       = data/x
                        newflat_name = f'{curdate}_n{band}{self.flattype}_bin{bin}.fits'
                        fits.writeto(newflat_name, nimage, header=newhdr, overwrite=True)

                        os.makedirs(f'{savedir}masterflat_{band}/', exist_ok=True)
                        shutil.copy(newflat_name, f'{savedir}masterflat_{band}/')
                        print('Normalised master flats are created.')
        iraf.imstat(images='*n*flat*.fits', fields='image,mean,midpt,mode,stddev,min,max', lsigma=3, usigma=3, Stdout=0)
        self.masterflat = glob.glob('*n*flat*.fits')
    '''
    def Apply(self) :
        """
        Apply bias and dark subtraction, and flat fielding to light frame images.

        This function processes light frame images by performing the following steps:
        1. Filter the summary table for 'Light Frame' images.
        2. Extract unique binning settings, exposure times, and filter bands.
        3. Perform bias subtraction on each light frame image.
        4. Perform dark subtraction on each light frame image.
        5. Apply flat fielding to each light frame image.
        6. Store the processed images in a specified directory.

        Attributes:
        - `summary`: A DataFrame containing metadata of images.
        - `curdate`: Current date used in naming the output files.
        - `savedir`: Directory where the processed light frame files will be saved.
        - `fdz_images`: List of fully processed light frame images created.

        Parameters:
        None

        Returns:
        None

        Raises:
        IOError: If there are issues with reading or writing image files.
        """
        sum       = self.summary
        curdate   = self.curdate
        savedir   = self.savedir

        sci           = sum[sum['IMAGETYP'] == 'Light Frame']
        scibinset     = list(set(sci['XBINNING'])) ; scibinset.sort()
        sciexpset     = list(set(sci['EXPTIME']))  ; sciexpset.sort()
        scibandset    = list(set(sci['FILTER']))   ; scibandset.sort()

        fdz_images = []
        for inim in sci['FILENAME'] :
            print(sci[sci['FILENAME'] == inim])
            imidx   = np.where(sci['FILENAME'] == inim)[0]
            bin     = sci['XBINNING'][imidx][0]
            exp     = sci['EXPTIME'][imidx][0]
            band    = sci['FILTER'][imidx][0]
            ccdtemp = sci['CCDTEMP0'][imidx][0]

            mblist = glob.glob(f'{self.bias_archive}2*_zero_bin{bin}.fits'); mblist.sort() 
            mdlist = glob.glob(f'{self.dark_archive}2*_dark{exp}_bin{bin}_{ccdtemp}.fits'); mdlist.sort() 
            mflist = glob.glob(f'{self.flat_archive}2*_n{band}skyflat_bin{bin}.fits'); mflist.sort() 


            if self.zeroproc :
                print('Zero subtraction...')
                #iraf.imarith(operand1=inim, op='-', operand2=f'./*_zero_bin{bin}.fits', result='z'+inim)
                iraf.imarith(operand1=inim, op='-', operand2=f"{self.curdir}/{mblist[-1].split('/')[-1]}", result='z'+inim)
            else :
                print('Zero subtraction is skipped...')
                os.system(f'cp {inim} z{inim}')

            if self.darkproc :
                print('Dark subtraction...')
                #iraf.imarith(operand1 = 'z'+inim, op = '-', operand2 = f'./*_dark{exp}_bin{bin}.fits', result = 'dz'+inim)
                iraf.imarith(operand1 = 'z'+inim, op = '-', operand2 = f"{self.curdir}/{mdlist[-1].split('/')[-1]}", result = 'dz'+inim)
            else :
                print('Dark subtraction is skipped...')
                os.system(f'cp z{inim} dz{inim}')

            if self.flatproc :
                print('Flat fielding...')
                dzinput = 'dz'+inim
                nflat   = glob.glob(f'2*n{band}{self.flattype}_bin{bin}.fits')[0]
                iraf.imarith(operand1=dzinput, op='/', operand2=nflat, result='f'+dzinput)      
            else:
                print('Flat fielding is skipped...!?')
                os.system(f'cp dz{inim} fdz{inim}')
            fdz_images.append(f'fdz{inim}')
        print('Flat fielding is finished. Check the images.') 
        self.fdz_images     = fdz_images
    '''
    def Apply(self):
        """
        Apply bias and dark subtraction, and flat fielding to light frame images.
        """
        sum     = self.summary
        curdate = self.curdate
        savedir = self.savedir

        sci = sum[sum['IMAGETYP'] == 'Light Frame']
        scibinset = list(set(sci['XBINNING']))
        scibinset.sort()
        sciexpset = list(set(sci['EXPTIME']))
        sciexpset.sort()
        scibandset = list(set(sci['FILTER']))
        scibandset.sort()

        fdz_images = []

        for inim in sci['FILENAME']:
            print(sci[sci['FILENAME'] == inim])
            imidx = np.where(sci['FILENAME'] == inim)[0]
            bin = sci['XBINNING'][imidx][0]
            exp = sci['EXPTIME'][imidx][0]
            band = sci['FILTER'][imidx][0]
            ccdtemp = sci['CCDTEMP0'][imidx][0]

            # Get master files
            bias_file = self.GetMaster('bias', bin=bin)
            dark_file = self.GetMaster('dark', bin=bin, exp=exp, temp=ccdtemp)
            flat_file = self.GetMaster('flat', bin=bin, filter_name=band)

            # Check for required files
            if not bias_file:
                raise IOError(f"No valid bias file found for bin: {bin}.")
            if not dark_file:
                raise IOError(f"No valid dark file found for bin: {bin}, exp: {exp}, temp: {ccdtemp}.")
            if not flat_file:
                raise IOError(f"No valid flat file found for bin: {bin}, band: {band}.")

            # Zero subtraction
            if self.zeroproc:
                print("Zero subtraction...")
                iraf.imarith(operand1=inim, op='-', operand2=bias_file, result='z' + inim)
            else:
                print("Zero subtraction is skipped...")
                os.system(f'cp {inim} z{inim}')

            # Dark subtraction
            if self.darkproc:
                print("Dark subtraction...")
                zinput = 'z' + inim
                iraf.imarith(operand1=zinput, op='-', operand2=dark_file, result='dz' + inim)
            else:
                print("Dark subtraction is skipped...")
                os.system(f'cp z{inim} dz{inim}')

            # Flat fielding
            if self.flatproc:
                print("Flat fielding...")
                dzinput = 'dz' + inim
                iraf.imarith(operand1=dzinput, op='/', operand2=flat_file, result='f' + dzinput)
            else:
                print("Flat fielding is skipped...")
                os.system(f'cp dz{inim} fdz{inim}')

            fdz_images.append(f'fdz{inim}')

        print("Flat fielding is finished. Check the images.")
        self.fdz_images = fdz_images

    def Astrometry(self, imlist_name = None) :
        """
        Perform astrometric calibration on processed light frame images using Astrometry.net.

        This function processes a list of fully processed light frame images (`fdz_images`)
        by performing the following steps:
        1. Read the configuration file for Source Extractor.
        2. Solve the World Coordinate System (WCS) for each image using Astrometry.net.
        3. Optionally, use RA and Dec from the image headers if available.
        4. Clean up temporary files created during the process.
        5. Store the successfully solved images in the `solved_images` attribute.

        Attributes:
        - `fdz_images`: List of fully processed light frame images.
        - `sexconfig`: Configuration file for Source Extractor.
        - `radius`: Search radius for astrometric calibration.
        - `solved_images`: List of successfully solved images with WCS information.

        Parameters:
        None

        Returns:
        None

        Raises:
        IOError: If there are issues with reading or writing image files.
        """

        if imlist_name == None:
            imlist      = self.fdz_images
        elif imlist_name != None:
            imlist   = glob.glob(imlist_name)

        sexconfig   = self.sexconfig
        radius      = self.radius
        scalelow    = self.scalelow
        scalehigh   = self.scalehigh

        os.system(f'cat {sexconfig}')
        print('Solving WCS using Astrometry.net...')
        for n, inim in enumerate(imlist) :
            if self.radec == False:
                com = f'solve-field {inim} --resort --cpulimit 120 --skip-solved --config {self.config} --use-source-extractor  --source-extractor-config {sexconfig} --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low {scalelow} --scale-high {scalehigh} --radius {radius} --no-remove-lines --uniformize 0 --no-plots  --new-fits a{inim} --overwrite --temp-dir .\n'
            else:
                hdr     = fits.getheader(inim)
                ra, dec = ':'.join(hdr['OBJCTRA'].split(' ')), ':'.join(hdr['OBJCTDEC'].split(' '))

                com = f'solve-field {inim} --resort --cpulimit 120 --skip-solved --config {self.config} --use-source-extractor  --source-extractor-config {sexconfig} --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low {scalelow} --scale-high {scalehigh} --radius {radius} --no-remove-lines --uniformize 0 --no-plots  --new-fits a{inim} --overwrite --temp-dir . --ra {ra} --dec {dec} --overwrite\n'

            print(com)
            print(f'{n} th of {len(imlist)}')
            os.system(com) 

        success = glob.glob('a*.fit') ; success.sort()
        print(f"From {len(imlist)} files, {len(success)} files are solved. ({len(success)/len(imlist)*100:.1f} %)")

        os.system('rm tmp*')
        os.system('rm *.wcs *.rdls *.corr *.xyls *.solved *.axy *.match')
        print('Astrometry process is complete.')
        self.solved_images = success

    def Check(self, imlist_name = None):
        """
        Verify and annotate solved images with target object information from a catalog.

        This function performs the following steps:
        1. Read observational specifications for the current camera.
        2. Read and process a catalog of target objects, converting RA and Dec to degrees.
        3. For each solved image, determine the image center's coordinates.
        4. Match the image center to the target catalog to identify the closest object.
        5. Check if the closest object is within the field of view.
        6. Annotate the image header with the matched object's name if within the field of view.

        Attributes:
        - `solved_images`: List of solved images with WCS information.

        Parameters:
        None

        Returns:
        None

        Raises:
        IOError: If there are issues with reading or writing image files.
        """
        if imlist_name == None:
            imlist = self.solved_images
        elif imlist_name != None:
            imlist   = glob.glob(imlist_name)

        all_catname   = self.alltarget
        all_cat       = ascii.read(all_catname) 
        ra, dec       = all_cat['ra'], all_cat['dec']
        radeg, decdeg = [], []
        for i in range(len(all_cat)) :
            c   = SkyCoord(str(ra[i])+' '+str(dec[i]), unit=(u.hourangle, u.deg))
            radeg.append(c.ra.deg)
            decdeg.append(c.dec.deg)
        all_cat['radeg']  = radeg
        all_cat['decdeg'] = decdeg
        coo_all = SkyCoord(radeg, decdeg, unit=(u.deg,u.deg))
        
        for inim in imlist :
            print(inim)
            data, hdr       = fits.getdata(inim, header=True)
            CRVAL1, CRVAL2  = wcscenter(inim)

            coo_target      = SkyCoord(CRVAL1, CRVAL2, unit=(u.deg, u.deg))
            indx, d2d, d3d  = coo_target.match_to_catalog_sky(coo_all)
            if d2d.arcmin > self.fov/2. :
                if ('TOI' or 'KOI' or 'TIC' or 'Qatar' or 'HAT' or 'KIC' or 'EPIC') in hdr['OBJECT'] :
                    print('Exoplanet observation...?')
                    hdr['PROGRAM'] = 'Exoplanet'
                elif ('SN' or 'AT' or 'ZTF' or 'ATLAS' or 'GRB' or 'Gaia') in hdr['OBJECT'] :
                    print('Transient observation...')
                    hdr['PROGRAM'] = 'Transient'
                else:
                    print('Unknown field.')
                    hdr['PROGRAM']  = ' '
                fits.writeto(inim, data, header=hdr, overwrite=True)
            elif d2d.arcmin <= self.fov/2. :
                obj = all_cat[indx]['obj']
                print('======================================')
                print(f'{obj} is matched in IMSNG field.')
                print(f'{round(d2d.arcmin[0],1)} arcmin apart')
                print('======================================')
                hdr['OBJECT']   = obj
                hdr['PROGRAM']  = 'IMSNG'
                fits.writeto(inim, data, header=hdr, overwrite=True)
        print('Headers are all checked.')

    def EditName(self, imlist_name = None) :
        """
        Rename solved images according to specified naming convention.

        This function renames solved images based on the following information extracted
        from the image header:
        - Exposure time (EXPTIME)
        - UT date and start time (UTDATE and UTSTART)
        - Filter used (FILTER)

        The new filename format is:
        'Calib-{ccd}-{OBJECT}-{UTDATE}-{UTSTART}-{FILTER}-{EXPTIME}.fits'

        Attributes:
        - `ccd`: The identifier for the camera or CCD used.
        - `solved_images`: List of solved images with WCS information.
        - `reduced_images`: List of renamed reduced images.

        Parameters:
        None

        Returns:
        None

        Raises:
        IOError: If there are issues with reading or writing image files.
        """
        if imlist_name == None:
            imlist = self.solved_images
            self.summary['CALNAME'] = self.summary['FILENAME'].copy()
            #self.summary['CALNAME'] = self.summary['CALNAME'].astype('str60') 
            newcol = Column(data=self.summary['FILENAME'], name='CALNAME', dtype=f'S60')
            self.summary.replace_column('CALNAME', newcol)
            for inim in imlist :
                idx     = np.where(self.summary['FILENAME'] == inim[4:])[0]
                hdr     = fits.getheader(inim)
                EXPTIME = int(hdr['exptime'])
                UTDATE  = hdr['date-obs'][0:10]
                UTSTART = hdr['date-obs'][11:19]
                OBJECT  = hdr['OBJECT']
                FILTER  = hdr['FILTER']

                newimage = f'Calib-{self.ccd}-{OBJECT}-{str(UTDATE[0:4])}{str(UTDATE[5:7])}{str(UTDATE[8:10])}-{str(UTSTART[0:2])}{str(UTSTART[3:5])}{str(UTSTART[6:8])}-{FILTER}-{str(EXPTIME)}.fits'
                #'''
                if self.summary['IMAGETYP'][idx].item() == 'Light Frame' :
                    self.summary['CALNAME'][idx] = newimage
                else:
                    self.summary['CALNAME'][idx] = inim
                #'''
                print(f'Copy {inim} to {newimage}...')
                os.system(f'/usr/bin/cp {inim} {newimage}')
            print(f"Basic preprocessing of {self.ccd} is finished.")

            print('Remove previous check images...')
            os.system(f'rm {self.curdir}/Cal*.ch*.fits')
            
        elif imlist_name != None:
            imlist   = glob.glob(imlist_name)
            for inim in imlist :
                #idx     = np.where(self.summary['FILENAME'] == inim[4:])[0]
                hdr     = fits.getheader(inim)
                EXPTIME = int(hdr['exptime'])
                UTDATE  = hdr['date-obs'][0:10]
                UTSTART = hdr['date-obs'][11:19]
                OBJECT  = hdr['OBJECT']
                FILTER  = hdr['FILTER']

                newimage = f'Calib-{self.ccd}-{OBJECT}-{str(UTDATE[0:4])}{str(UTDATE[5:7])}{str(UTDATE[8:10])}-{str(UTSTART[0:2])}{str(UTSTART[3:5])}{str(UTSTART[6:8])}-{FILTER}-{str(EXPTIME)}.fits'

                print(f'Copy {inim} to {newimage}...')
                os.system(f'/usr/bin/cp {inim} {newimage}')
            print(f"Basic preprocessing of {self.ccd} is finished.")  
        self.reduced_images = glob.glob('Cal*.fits') ; self.reduced_images.sort()

    def FluxScaling(self, imlist_name = None, zp0=25, zpkey='ZP'):

        if imlist_name == None:
            imlist      = self.reduced_images 
        elif imlist_name != None:
            imlist   = glob.glob(imlist_name); imlist.sort()
            
        for inim in imlist :
            print(f'Flux scaling for {inim}...')
            data, hdr = fits.getdata(inim, header=True)
            if 'ZP' in hdr : 
                zp = float(hdr[zpkey])
            else:
                print(f'{inim}: no ZP in header. PASS.')
                pass
	
            factor_dum = 10.**(0.4*(zp0 - zp))
            print(f'Scale factor = {factor_dum:.2f}')

            newdata_name    = 'f'+inim[:-5]+'.fits'
            newdata         = data * factor_dum
            hdr['FLXSCALE'] = factor_dum 
            print(f'Writing new data as {newdata_name}...')
            fits.writeto(newdata_name, newdata, hdr, overwrite=True)

        self.fscaled_images = glob.glob('fCal*.fits'); self.fscaled_images.sort()
        print('Done.')

    def Stack(self, imlist_name = None):
        """
        Stack calibrated light frame images by object, filter, and binning.

        This function performs the following steps:
        1. Identify unique objects, filters, and binning settings from the list of calibrated light frames.
        2. For each combination of object, filter, and binning:
            a. Select the relevant images.
            b. Group images by epoch using the `SameEpoch` function.
            c. Perform WCS remapping on the images using the `run_wcsremap_epoch` function.
            d. Combine the remapped images using the `imcombine_set` function.
        3. Handle cases where no images are found for a particular combination.
        4. Print status messages throughout the process.

        Attributes:
        - `sep`: Separator used in the `SameEpoch` function.
        - `wrapper`: Wrapper used for WCS registration in the `run_wcsremap_epoch` function.
        - `reject`: Rejection method used in the `imcombine_set` function.

        Parameters:
        None

        Returns:
        None

        Raises:
        IOError: If there are issues with reading or writing image files.
        """
        #imlist_name = 'Cal*.fits'
        #imlist_name = 'a*.fit'
      
        if imlist_name == None:
            imlist      = self.fscaled_images
            objset  =  list(set(self.summary['OBJECT'][self.summary['IMAGETYP'] == 'Light Frame']))
            bandset =  list(set(self.summary['FILTER'][self.summary['IMAGETYP'] == 'Light Frame']))
            binset  =  list(set(self.summary['XBINNING'][self.summary['IMAGETYP'] == 'Light Frame']))
           

        elif imlist_name != None:
            #imlist   = glob.glob(imlist_name); imlist.sort()

            objset  = list(set(iraf.hselect(imlist_name, fields="$OBJECT", expr=f'IMAGETYP = "Light Frame"', Stdout = 1))) ; objset.sort()
            bandset = list(set(iraf.hselect(imlist_name, fields="$FILTER", expr=f'IMAGETYP = "Light Frame"', Stdout = 1))) ; bandset.sort()
            binset  = list(set(iraf.hselect(imlist_name, fields="$XBINNING", expr=f'IMAGETYP = "Light Frame"', Stdout = 1))) ; binset.sort()

        for bin in binset :
            for obj in objset :
                for band in bandset:
                    if imlist_name == None:
                        calframe = list(self.summary['CALNAME'][(self.summary['FILTER'] == band) & (self.summary['IMAGETYP'] == 'Light Frame')  & (self.summary['XBINNING'] == bin) & (self.summary['OBJECT'] == obj)])
                    elif imlist_name != None:
                        calframe = iraf.hselect(imlist_name, fields="$I", expr=f'OBJECT == "{obj}" && IMAGETYP == "Light Frame" && FILTER == "{band}" && XBINNING == {bin}', Stdout = 1)
                    print(calframe)
                    if len(calframe) != 0:
                        try:
                            image_list = SameEpoch(f"fCal*{obj}*-{band}-*.fits", sep=self.sep)
                            print(image_list)

                            imlist = list(image_list[0])
                            refim  = imlist[0]
                            
                            run_wcsremap_epoch(f"fCal*{obj}*-{band}-*.fits", sep=self.sep, wcsregister=self.wrapper)
                            #imcombine_set(f"Remap*-{obj}-*-{band}-*.fits", reject=self.reject)
                            imcombine_epoch(f"Remap*-{obj}-*-{band}-*.fits", reject=self.reject, sep=self.sep)
                        except:
                            print(f'No fCal*{obj}*-{band}-*.fits or error in wcsremap.')
                            continue
        print('Finished.')
        self.combined_images = glob.glob('Cal*com.fits'); self.combined_images.sort()

    def Dophot(self, imlist_name = None, ref='PS1', magup=12, maglow=17):
        """
        Single image와 combined image 모두 측광해둘 예정
        """

        if imlist_name == None:
            imlist      = self.reduced_images + self.combined_images 
        elif imlist_name != None:
            imlist   = glob.glob(imlist_name); imlist.sort()

        if not os.path.exists(f'./photretry') :
            os.mkdir(f'./photretry')

        for inim in imlist :
            print(f'Do photomtry for {inim}')
            hdr     = fits.getheader(inim)
            obj     = hdr['OBJECT']
            ccd     = self.ccd
            band    = hdr['FILTER']
            program = hdr.get('PROGRAM', '')

            if band in ['B', 'V', 'R', 'I', 'g','r','i','z'] :

                if ref == None:
                    try:
                        phot(inim, target_catalog='', band=band, path='/home/lim9/miniconda3/lib/python3.9/site-packages/lgpytars/photconf/', savecat=True, subprefix='hdCalib', ref='PS1', mykey='APER_1', sub=False, minarea=3, det_thresh=3, deb_nthresh=32, deb_mincont=0.01, backsize=128, backfiltersize = 5, backphoto_type='LOCAL',  ratio=1, onlyzp=True, snrcut=0.1, magup=12, maglow=17)

                    except:
                        try:
                            phot(inim, target_catalog='', band=band, path='/home/lim9/miniconda3/lib/python3.9/site-packages/lgpytars/photconf/', savecat=True, subprefix='hdCalib', ref='APASS', mykey='APER_1', sub=False, minarea=3, det_thresh=3, deb_nthresh=32, deb_mincont=0.01, backsize=128, backfiltersize = 5, backphoto_type='LOCAL', ratio=1, onlyzp=True, snrcut=0.1, magup=11, maglow=15)
                        except:
                            print('No stars matched with APASS or PS1. PLEASE DO MANUALLY.')
                        
                            os.system(f'mv {inim} ./photretry/')
                            os.system(f'mv {inim[:-5]}.chaper.fits ./photretry/')
                            os.system(f'mv {inim[:-5]}.chbkg.fits ./photretry/')
                            os.system(f'mv {inim[:-5]}.chseg.fits ./photretry/')
                            os.system(f'mv {inim[:-5]}.merge.cat ./photretry/')
                            print(f'{inim} is moved to photretry directory.')

                            self.reduced_images.remove(inim)
                else:
                    phot(inim, target_catalog='', band=band, path='/home/lim9/miniconda3/lib/python3.9/site-packages/lgpytars/photconf/', savecat=True, subprefix='hdCalib', ref=ref, mykey='APER_1', sub=False, minarea=3, det_thresh=3, deb_nthresh=32, deb_mincont=0.01, backsize=128, backfiltersize = 5, backphoto_type='LOCAL',  ratio=1, onlyzp=True, snrcut=0.1, magup=magup, maglow=maglow)
                        
            else:
                print(f'{band}-band is not supported. TBD.')
        os.system('ls ./photretry/*.fits')
        print('Finished.')  

    def Subtract(self, imlist_name='Cal*.com.fits', template_dir=None):
        """
        Perform image subtraction using HOTPANTS.
        
        This function identifies the appropriate template image for each science image,
        then uses HOTPANTS to perform subtraction.

        Parameters:
        imlist_name (str): Pattern to match science images.
        template_dir (str): Directory containing the template images.

        Returns:
        None
        """


        if imlist_name == None:
            imlist      = self.combined_images 
        elif imlist_name != None:
            imlist   = glob.glob(imlist_name); imlist.sort()

        if template_dir == None :
            template_dir = f'{self.template_dir}/{self.ccd}'
            #self.template_dir = template_dir
        elif template_dir != None:
            pass
            print(f'Template path is : {template_dir}')

        # Iterate through science images
        no_temp = []
        for inim in imlist:
            hdr  = fits.getheader(inim)
            obj  = hdr['OBJECT']
            band = hdr['FILTER']

            if hdr['program'] != 'IMSNG':
                print(f'{inim} is not IMSNG program.')
                pass
            
            # Match template image
            template_pattern = f"{template_dir}/Calib-PNUO_C361K-{obj}-*-{band}-*.fits"
            templates = glob.glob(template_pattern)
            
            if not templates:
                print(f"No matching template found for {inim}. Skipping...")
                no_temp.append(inim)
                continue
            
            # Select the best matching template (e.g., closest date if multiple templates exist)
            template = templates[0]  # Simplified; can be expanded for best match logic
            print(f"Using template: {template} for {inim}")

            # Registration with template
            print(f"Aligning {template} into {inim}...")
            rtemp = run_wcsremap_wrap(inim, template, large=True, wcsregister=True)
            
            # Perform subtraction with HOTPANTS
            hotpants.runhotpants(imlist_name=inim, refinst = self.ccd, scifwhmkey='FWHM', reffwhmkey='FWHM', refim='', same=True, part0='Remap')
        
        print("Subtraction complete for all matching images.")
        print('Images with no templates')
        print(no_temp)

        self.sub_images = glob.glob('hd*com.fits'); self.sub_images.sort()
        #self.conv_images = glob.glob('hc*com.fits'); self.conv_images

    def Archive(self, imlist_name = None):
        """
        Archive processed images into designated storage based on the CCD and program type.

        Parameters:
        - imlist_name (str, optional): If provided, only images matching this pattern will be archived.

        The method:
        - Determines the archive path based on self.ccd and program type (IMSNG or local).
        - Creates necessary directories if they don't exist.
        - Moves the processed images to the appropriate archive.
        """
        
        #imsng_archive = '/mnt/dataset/obsdata/IMSNG' 
        #pnuo_archive  = '/mnt/dataset/obsdata/PNUO'
        archive_base = {
            "PNUO": "/mnt/dataset/obsdata/PNUO",
            "MAAO": "/mnt/dataset/obsdata/MAAO",
            #"OTHER": "/mnt/dataset/obsdata/OTHER"  # Default directory for unknown CCDs
        }

        imsng_archive     = "/mnt/dataset/obsdata/IMSNG"
        transient_archive = "/mnt/dataset/obsdata/Transient"
        exoplanet_archive = "/mnt/dataset/obsdata/Exoplanet"

        if imlist_name is None:
            imlist = self.reduced_images + self.combined_images
            if hasattr(self, 'sub_images') and self.sub_images is not None:
                imlist += self.sub_images 
            self.imlist = imlist
        elif imlist_name is not None:
            imlist = glob.glob(imlist_name); imlist.sort()

        for inim in imlist:
            hdr     = fits.getheader(inim)
            obj     = hdr['OBJECT']

            #ccd     = self.ccd
            band    = hdr['FILTER']
            program = hdr.get('PROGRAM', '')
            #part    = inim[:-5]
            # others  = glob.glob(part+'*')
            #"""
            if program == 'IMSNG':
                print(f'Move {inim} to {imsng_archive}')
                archive_path = os.path.join(imsng_archive, obj, self.ccd, band)
            elif program == 'Transient':
                print(f'Move {inim} to {transient_archive}')
                archive_path = os.path.join(transient_archive, obj, self.ccd, band)     
            elif program == 'Exoplanet' :
                print(f'Move {inim} to {exoplanet_archive}')
                archive_path = os.path.join(exoplanet_archive, obj, self.ccd, band)                             
            else:
                if self.ccd == 'MAAO_STX16803':
                    print(f"Move {inim} to {archive_base['MAAO']}")
                    archive_path = os.path.join(archive_base['MAAO'], obj, self.ccd, band)
                elif self.ccd == 'PNUO_C361K': 
                    print(f"Move {inim} to {archive_base['PNUO']}")
                    archive_path = os.path.join(archive_base['PNUO'], obj, self.ccd, band)
            #"""

            os.makedirs(archive_path, exist_ok=True)
            dest_file_path = os.path.join(archive_path, os.path.basename(inim))
            os.system(f'/usr/bin/cp {inim} {dest_file_path}')

            os.system(f'chmod -R 777 {archive_path}')
        print('Finished.')

    #'''
#if __name__ == "__main__":
def main(ccd='PNUO_C361K'):
    # Create an instance of the Red class
    
    red = Red(imlist_name = '*.fit',
                sumfile     = 'file_summary.txt',
                ccd         = ccd,
                zeroproc    = True, 
                darkproc    = True,
                flatproc    = True,
                flattype    = 'skyflat'
                )
    

    # Execute the methods in the desired sequence
    if os.path.exists(red.sumfile):
        print(f'{red.sumfile} already exists.')
        red.summary = ascii.read(red.sumfile)
        red.summary.sort('MJD')       
    else:     
        red.FileSum()

    imtype = list(set(red.summary['IMAGETYP']))
    if 'Bias Frame' in imtype:
        red.GenBias()
    else:
        print('Master bias exists already or no bias in archive.')

    if 'Dark Frame' in imtype:
        red.GenDark()
    else:
        print('Master dark exists already or no dark in archive.')

    if 'Flat Field' in imtype:
        red.GenFlat()
    else:
        print('Master flat exists already or no flat in archive.')

    if 'Light Frame' in imtype:
        red.Apply()
    else:
        print('Science images exist already or do not exist.')
    
    if hasattr(red, 'fdz_images'):
        red.Astrometry()
        red.Check()
        red.EditName()
    else:
        print('No reduced images here.')

    if hasattr(red, 'reduced_images'):
        os.system(f'rm {red.curdir}/Cal*.ch*.fits')
        red.Dophot(imlist_name = 'Cal*0.fits')
    else:
        print('No single images to do photometry.')

    if hasattr(red, 'reduced_images') :
        red.FluxScaling(imlist_name=None, zp0=25, zpkey='ZP')
    else:
        print('No images to do flux scaling.')

    if hasattr(red, 'fscaled_images') :
        red.Stack(imlist_name=None)
    else:
        print('No images to stack.')

    if hasattr(red, 'combined_images') :
        red.Dophot(imlist_name = 'Cal*0.com.fits')
    else:
        print('No stacked images to do photometry.')

    if hasattr(red, 'combined_images') :
        red.Subtract(imlist_name = None, template_dir=None)

    red.Archive()

    print('gethead Cal*0.com.fits fwhm zp zper lim5 lim3')
    os.system('gethead Cal*0.com.fits fwhm zp zper lim5 lim3')
    print(f'{red.ccd} data reduction ({red.curdate}) is complete.')
#'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LGPYRED data reduction pipeline")
    parser.add_argument("--all", action="store_true", help="Run the entire pipeline")
    parser.add_argument("--filesum", action="store_true", help="Run FileSum")
    parser.add_argument("--genbias", action="store_true", help="Run GenBias")
    parser.add_argument("--gendark", action="store_true", help="Run GenDark")
    parser.add_argument("--genflat", action="store_true", help="Run GenFlat")
    parser.add_argument("--apply", action="store_true", help="Run Apply")
    parser.add_argument("--astrometry", action="store_true", help="Run Astrometry")
    parser.add_argument("--check", action="store_true", help="Run Check")
    parser.add_argument("--editname", action="store_true", help="Run EditName")
    parser.add_argument("--dophot", action="store_true", help="Run Dophot")
    parser.add_argument("--fluxscaling", action="store_true", help="Run FluxScaling")
    parser.add_argument("--stack", action="store_true", help="Run Stack")
    parser.add_argument("--subtract", action="store_true", help="Run Subtract")
    parser.add_argument("--archive", action="store_true", help="Run Archive")
    parser.add_argument("--imlist", type=str, default=None, help="Specify image list pattern (e.g., 'Cal*com.fits')")
    parser.add_argument("--template_dir", type=str, default="/mnt/dataset/obsdata/IMSNG/template_20250213/", help="Directory for template images")

    args = parser.parse_args()

    # Create an instance of the Red class
    red = Red(imlist_name=args.imlist if args.imlist else '*.fit',
              sumfile='file_summary.txt',
              ccd='PNUO_C361K',
              zeroproc=True,
              darkproc=True,
              flatproc=True,
              flattype='skyflat')

    # Execute commands based on arguments
    if args.all:
        args.filesum = args.genbias = args.gendark = args.genflat = args.apply = True
        args.astrometry = args.check = args.editname = args.dophot = True
        args.fluxscaling = args.stack = args.subtract = args.archive = True

    if args.filesum:
        if os.path.exists(red.sumfile):
            print(f'{red.sumfile} already exists.')
            red.summary = ascii.read(red.sumfile)
            red.summary.sort('MJD')
        else:
            red.FileSum()

    if args.genbias:
        imtype = list(set(red.summary['IMAGETYP']))
        if 'Bias Frame' in imtype:
            red.GenBias()
        else:
            print('Master bias exists already or no bias in archive.')

    if args.gendark:
        imtype = list(set(red.summary['IMAGETYP']))
        if 'Dark Frame' in imtype:
            red.GenDark()
        else:
            print('Master dark exists already or no dark in archive.')

    if args.genflat:
        imtype = list(set(red.summary['IMAGETYP']))
        if 'Flat Field' in imtype:
            red.GenFlat()
        else:
            print('Master flat exists already or no flat in archive.')

    if args.apply:
        imtype = list(set(red.summary['IMAGETYP']))
        if 'Light Frame' in imtype:
            red.Apply()
        else:
            print('Science images exist already or do not exist.')

    if args.astrometry and hasattr(red, 'fdz_images'):
        red.Astrometry(imlist_name=args.imlist)

    if args.check and hasattr(red, 'fdz_images'):
        red.Check(imlist_name=args.imlist)

    if args.editname and hasattr(red, 'fdz_images'):
        red.EditName(imlist_name=args.imlist)

    if args.dophot and hasattr(red, 'reduced_images'):
        red.Dophot(imlist_name=args.imlist)

    if args.fluxscaling and hasattr(red, 'reduced_images'):
        red.FluxScaling(imlist_name=args.imlist, zp0=25, zpkey='ZP')

    if args.stack and hasattr(red, 'fscaled_images'):
        red.Stack(imlist_name=args.imlist)

    if args.dophot and hasattr(red, 'combined_images'):
        red.Dophot(imlist_name = None)
    else:
        print('No stacked images to do photometry.')

    if args.subtract and hasattr(red, 'combined_images'):
        red.Subtract(imlist_name=args.imlist, template_dir=None)

    if args.archive:
        red.Archive(imlist_name=args.imlist)

    print(f'{red.ccd} data reduction ({red.curdate}) is complete.')
