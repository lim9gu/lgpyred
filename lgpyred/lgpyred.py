### lgpyred image REDuction pipeline

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

from importlib.resources import files
import importlib.resources as pkg_resources

import lgpyred.data
from lgpyred.SameEpoch import SameEpoch
from lgpyred.imcombine import imcombine_set, imcombine_epoch
from lgpyred.reduction.hdrcheck import wcscenter
from lgpyred.reduction.wcsremap import run_wcsremap_epoch, run_wcsremap_wrap
from lgpyred.lgphot import phot
from lgpyred.reduction import hotpants

class Red :
    """
    Red: A CCD image reduction pipeline for astronomical observations.

    This class performs bias, dark, flat corrections, WCS alignment, image stacking,
    image subtraction using HOTPANTS, and optional photometry tasks. It supports 
    multiple observatories and CCD types defined in an instrument specification table.

    It also supports archiving reduced data to user-defined directories, and references
    external configuration files such as astrometry.net and SExtractor config.

    Note:
    - IRAF must be installed and accessible.
    - Astrometry.net and HOTPANTS must be available in the system PATH.
    - Configuration and reference files (obs_spec.txt, alltarget.dat, sex config) are loaded from lgpyred package.
    """
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
        Initialize the Red class for CCD image reduction and analysis.

        This class is designed to automatically perform CCD image calibration (bias, dark, flat),
        astrometric solving, image stacking, subtraction, and archiving. Instrument specifications 
        are loaded based on the specified CCD from packaged data files.

        Parameters
        ----------
        imlist_name : str, optional
            Pattern to match image files to reduce (default: '*.fit').

        sumfile : str, optional
            Output file to store image summary information (default: 'file_summary.txt').

        ccd : str, optional
            CCD name used to fetch instrument-specific configuration (default: 'PNUO_C361K').

        zeroproc : bool, optional
            Whether to apply bias correction (default: True).

        darkproc : bool, optional
            Whether to apply dark current correction (default: True).

        flatproc : bool, optional
            Whether to apply flat field correction (default: True).

        flattype : str, optional
            Type of flat field to use: 'skyflat' or 'Domeflat' (default: 'skyflat').

        sep : float, optional
            Time interval (in minutes) used for automatic image alignment (default: 60.0).

        archive_paths : dict, optional
            Custom dictionary to define archive destination directories. If None, default paths are used.

        Notes
        -----
        - Instrument parameters are loaded from the packaged `obs_spec.txt` and `ccddb/`.
        - All necessary configuration and support files (e.g., alltarget.dat, sextractor config) are 
          accessed via the `lgpyred.data`, `lgpyred.astrom_config`, and `lgpyred.photconf` modules.
        - If the `ASTROMETRY_CFG` environment variable is set, it overrides the default astrometry config path.
        - If the HOTPANTS template directory does not exist, subtraction will be skipped.
        - Archiving is automatically handled based on the `PROGRAM` header, CCD name, and filter band.

        Raises
        ------
        ValueError
            If `flattype` is not one of 'skyflat' or 'Domeflat'.
        """
        # File setting
        self.curdir      = os.getcwd()

        iraf.chdir(self.curdir)
        self.irafdb = str(files(lgpyred.data).joinpath("ccddb"))+'/'
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
        obscat           = ascii.read(str(files(lgpyred.data).joinpath('obs_spec.txt')))
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

        # Reduction setting
        self.zeroproc    = zeroproc
        self.darkproc    = darkproc
        self.flatproc    = flatproc

        self.flat_zeroproc = True
        self.flat_darkproc = True

        if flattype not in ['skyflat', 'Domeflat'] :
            raise ValueError(f"Invalid flattype: {flattype}. Must be 'skyflat' or 'Domeflat'.")
        self.flattype    = flattype

        self.bias_archive   = f'{self.savedir}masterbias/'
        self.dark_archive   = f'{self.savedir}masterdark/'
        self.flat_archive   = f'{self.savedir}masterflat_*/'

        # Astrometry setting (Astrometry.net)
        self.config    = os.environ.get('ASTROMETRY_CFG', '/etc/astrometry.cfg')

        self.scalelow  = self.pixscale * 0.9 #0.18  # lower limit for the pixel scale covering 'No binning' to '2x2 binning')
        self.scalehigh = self.pixscale * 2 * 1.1 #0.46  # upper limit for the pixel scale covering 'No binning' to '2x2 binning')
        self.radec     = False # if True, astronometry.net search for wcs using radec in fits image header 

        self.sexconfig = str(files('lgpyred.astrom_config'))

        self.radius    = 0.7 # degree of querying radius of stars in index files

        # Header check
        self.alltarget = str(files('lgpyred.data').joinpath('alltarget.dat'))

        # Image alignment setting (wcsremap)
        self.sep         = sep   # mins time interval
        self.wrapper     = True # use wcsremap wrapper?

        # Image stacking setting (IRAF.imcombine)
        self.reject = 'none' 

        # Photometry
        self.photpath = str(files('lgpyred.photconf'))

        # Image subtraction setting (HOTPANTS)
            # If you want to perform image subtraction, you need to prepare template images beforehand.
            # Please place your template images in the directory below, or change it to your own path.
            # If the directory does not exist, template_dir will be set to None and subtraction will be skipped.
        self.template_dir = '/mnt/dataset/obsdata/IMSNG/template_20250213'

        # Archive base directories (customize if needed)
        self.archive_paths = {
            'IMSNG':     '/mnt/dataset/obsdata/IMSNG',
            'Transient': '/mnt/dataset/obsdata/Transient',
            'Exoplanet': '/mnt/dataset/obsdata/Exoplanet',
            'PNUO':      '/mnt/dataset/obsdata/PNUO',
            'MAAO':      '/mnt/dataset/obsdata/MAAO',
            #'OTHER':   '/mnt/dataset/obsdata/OTHER'  # Add as needed
        }

        # File transfer setting
        self.server = ''

        self.banner()

        print(f'Reduction for {self.curdir} is now ready.')

    #@staticmethod
    def banner(self):
        """
        Print an ASCII banner displaying the pipeline name, version, and CCD.

        This method is called automatically during initialization to visually indicate
        that the LGPYRED pipeline is ready, along with the selected CCD and release version.

        The banner includes:
        - The name of the pipeline ("lgpyred imaging REDuction")
        - The current version
        - The CCD identifier (e.g., 'PNUO_C361K')

        No parameters or return values.
        """
        print(rf"""
        ██╗      ██████╗ ██████╗ ██╗   ██╗██████╗ ███████╗██████╗ 
        ██║     ██╔════╝ ██╔══██╗╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
        ██║     ██║  ███╗██████╔╝ ╚████╔╝ ██████╔╝█████╗  ██║  ██║
        ██║     ██║   ██║██╔═══╝   ╚██╔╝  ██╔══██╗██╔══╝  ██║  ██║
        ███████╗╚██████╔╝██║        ██║   ██║  ██║███████╗██████╔╝
        ╚══════╝ ╚═════╝ ╚═╝        ╚═╝   ╚═╝  ╚═╝╚══════╝╚═════╝ 
        """)
        print(rf"""
        lgpyred imaging REDuction pipeline for {self.ccd} 
        Release version: v1.0.0 (2025-07-07) 
        Author: G. Lim (lim9gu@gmail.com)
        """)
        # banner from https://patorjk.com/software/taag/#p=display&v=0&f=ANSI%20Shadow&t=lgpyred%0A-PNUO

    def GetMaster(self, filetype, bin=None, exp=None, temp=None, filter_name=None):
        """
        Locate the closest master calibration file (bias, dark, or flat) for the given observation date.

        Parameters:
        ----------
        imtype : str
            Type of calibration file to retrieve ('zero', 'dark', or 'flat').
        obsdate : str, optional
            Observation date in 'YYYYMMDD' format. If None, uses the current folder name.

        Returns:
        -------
        master_file : str
            Full path to the closest matching master calibration file.

        Notes:
        -----
        - Searches for calibration files within the archive directory defined in `self.bias_archive`,
        `self.dark_archive`, or `self.flat_archive`, depending on `imtype`.
        - Selects the calibration file with the closest date that is not after the observation date.
        - If no suitable file is found, returns an empty string.
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
        Generate a summary table for the given list of FITS images.

        Parameters:
        ----------
        imlist : list of str, optional
            List of image filenames. If None, matches files using `self.imlist_name`.

        Actions:
        -------
        - Extracts basic header information such as OBJECT, FILTER, EXPTIME, CCD-TEMP, etc.
        - Saves the summary as an ASCII table to `self.sumfile`.
        - Stores the table in `self.filetable`.

        Returns:
        -------
        None
        """
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
                    shutil.copy2(closest_bias, './')

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
                    shutil.copy2(closest_dark, './')
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
                    shutil.copy2(closest_flat, './')
                else:
                    print("No flat available.")
        print('Inspection finished.')

        self.summary = cat
        self.summary.sort('MJD')
    
    def GenBias(self):
        """
        Generate a master bias frame by combining zero-exposure images.

        Actions:
        --------
        - Identifies all images with IMAGETYP='zero' from the file summary.
        - Creates a master bias frame using IRAF's `zerocombine` task.
        - Saves the master bias frame in the archive directory defined by `self.bias_archive`.
        - The filename follows the convention: 'Bias_<curdate>.fits'.

        Returns:
        --------
        masterbias : str
            The filename of the generated master bias image.
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

    def GenDark(self) :
        """
        Generate a master dark frame by combining dark-exposure images.

        Actions:
        --------
        - Identifies all images with IMAGETYP='dark' from the file summary.
        - Subtracts the master bias (generated via `GenBias`) from each dark frame.
        - Creates a master dark using IRAF's `darkcombine` task.
        - Saves the master dark frame in the archive directory defined by `self.dark_archive`.
        - The filename follows the convention: 'Dark_<curdate>.fits'.

        Returns:
        --------
        masterdark : str
            The filename of the generated master dark image.
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
        Generate a master flat frame by combining skyflats or domeflats.

        Actions:
        --------
        - Identifies all images with IMAGETYP='flat' from the file summary.
        - Subtracts master bias and/or dark if `self.flat_zeroproc` or `self.flat_darkproc` is True.
        - Normalizes individual flat images.
        - Combines them into a master flat using IRAF's `flatcombine` task.
        - Saves the master flat frame in the archive directory defined by `self.flat_archive`.
        - The filename follows the convention: 'Flat_<filter>_<curdate>.fits'.

        Returns:
        --------
        masterflat : str
            The filename of the generated master flat image.
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
                iraf.imarith(operand1=im, op='-', operand2=bias_file, result='z' + im)
            else :
                print('Zero subtraction is skipped...')
                shutil.copy2(im, f'z{im}')

            if self.flat_darkproc :
                print('Dark subtraction...')
                zinput = 'z'+im
                iraf.imarith(operand1=zinput, op='-', operand2=dark_file, result='d' + zinput)

            else :
                print('Dark subtraction is skipped...')
                shutil.copy2(f'z{im}', f'dz{im}')

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

                        print(f'output : {output_name}')
                        iraf.imstat(images='@'+input_name,fields='image,mean,midpt,mode,stddev,min,max', lsigma=3, usigma=3)
                        iraf.flatcombine(input='@'+input_name, output=output_name, combine='median', reject='avsigclip', process='no', scale='mode', ccdtype='', lsigma='3.', hsigma='3.')

                        print(f'{output_name} is created. Normalizing...')
                        data, newhdr = fits.getdata(output_name, header=True)
                        x            = np.median(data)
                        nimage       = data/x
                        newflat_name = f'{curdate}_n{band}{self.flattype}_bin{bin}.fits'
                        fits.writeto(newflat_name, nimage, header=newhdr, overwrite=True)

                        os.makedirs(f'{savedir}masterflat_{band}/', exist_ok=True)
                        shutil.copy(newflat_name, f'{savedir}masterflat_{band}/')
                        print('Normalised master flats are created.')
        iraf.imstat(images='*n*flat*.fits', fields='image,mean,midpt,mode,stddev,min,max', lsigma=3, usigma=3, Stdout=0)
        self.masterflat = glob.glob('*n*flat*.fits')

    def Apply(self):
        """
        Apply calibration frames (bias, dark, flat) to science images.

        Actions:
        --------
        - Identifies all science images with IMAGETYP='object' from the file summary.
        - Performs overscan and trimming correction if applicable.
        - Applies master bias subtraction if `self.zeroproc` is True.
        - Applies master dark subtraction if `self.darkproc` is True.
        - Applies flat field correction if `self.flatproc` is True.
        - Saves reduced science frames with suffix `_red.fits`.

        Output:
        -------
        - Updates self.reduced_images with the list of processed science images.
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
                shutil.copy2(inim, f'z{inim}')

            # Dark subtraction
            if self.darkproc:
                print("Dark subtraction...")
                zinput = 'z' + inim
                iraf.imarith(operand1=zinput, op='-', operand2=dark_file, result='dz' + inim)
            else:
                print("Dark subtraction is skipped...")
                shutil.copy2(f'z{inim}', f'dz{inim}')

            # Flat fielding
            if self.flatproc:
                print("Flat fielding...")
                dzinput = 'dz' + inim
                iraf.imarith(operand1=dzinput, op='/', operand2=flat_file, result='f' + dzinput)
            else:
                print("Flat fielding is skipped...")
                shutil.copy2(f'dz{inim}', f'fdz{inim}')

            fdz_images.append(f'fdz{inim}')

        print("Flat fielding is finished. Check the images.")
        self.fdz_images = fdz_images

    def Astrometry(self, imlist_name = None) :
        """
        Solve astrometry for reduced science images using Astrometry.net.

        Actions:
        --------
        - Loops through reduced science images in `self.reduced_images`.
        - Calls Astrometry.net via `solve-field` command-line interface.
        - Uses custom configuration: pixel scale bounds (`scalelow`, `scalehigh`), 
        radius, and `sexconfig` file for source extraction.
        - Produces WCS-calibrated FITS files with `.wcs.fits` extension.

        Output:
        -------
        - Updates image headers with WCS solution if successful.
        - Skips images that already contain valid WCS headers.

        Notes:
        ------
        - Requires astrometry.net and index files to be properly installed on the system.
        - Configuration file can be set via the ASTROMETRY_CFG environment variable.
        """

        if imlist_name == None:
            imlist      = self.fdz_images
        elif imlist_name != None:
            imlist   = glob.glob(imlist_name)

        self.sexconfig = str(files('lgpyred.astrom_config'))

        sexconfig   = self.sexconfig
        radius      = self.radius
        scalelow    = self.scalelow
        scalehigh   = self.scalehigh

        sexparam    = str(files('lgpyred.astrom_config').joinpath('astrometry.net.param'))
        sexconv     = str(files('lgpyred.astrom_config').joinpath('astrometry.net.conv'))
        sexnnw      = str(files('lgpyred.astrom_config').joinpath('astrometry.net.nnw'))

        os.system(f'cat {sexconfig}')
        print('Solving WCS using Astrometry.net...')
        for n, inim in enumerate(imlist) :
            if self.radec == False:
                com = f"solve-field {inim} --resort --cpulimit 300 --skip-solved --config {self.config} --use-source-extractor --source-extractor-path 'source-extractor -c {self.sexconfig} -PARAMETERS_NAME {sexparam} -FILTER_NAME {sexconv} -STARNNW_NAME {sexnnw}' --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low {scalelow} --scale-high {scalehigh} --radius {radius} --no-remove-lines --uniformize 0 --no-plots  --new-fits a{inim} --overwrite --temp-dir .\n"
            else:
                hdr     = fits.getheader(inim)
                ra, dec = ':'.join(hdr['OBJCTRA'].split(' ')), ':'.join(hdr['OBJCTDEC'].split(' '))

                com = f"solve-field {inim} --resort --cpulimit 300 --skip-solved --config {self.config} --use-source-extractor --source-extractor-path 'source-extractor -c {self.sexconfig} -PARAMETERS_NAME {sexparam} -FILTER_NAME {sexconv} -STARNNW_NAME {sexnnw}' --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_AUTO --sort-ascending --scale-unit arcsecperpix --scale-low {scalelow} --scale-high {scalehigh} --radius {radius} --no-remove-lines --uniformize 0 --no-plots  --new-fits a{inim} --overwrite --temp-dir . --ra {ra} --dec {dec} --overwrite\n"

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
        Perform basic header checks and verification for science images.

        Actions:
        --------
        - Checks all science images for required FITS header keywords:
        OBJECT, DATE-OBS, FILTER, EXPTIME, CCD-TEMP, etc.
        - Replaces missing or invalid values with defaults (e.g., CCD-TEMP with self.temp2replace).
        - Compares target name in header against a reference list (`self.alltarget`)
        to validate naming and consistency.

        Output:
        -------
        - Logs or prints warnings for missing or incorrect header values.
        - Ensures all relevant headers are present before further reduction steps.
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
        Rename FITS files to follow a standardized naming convention.

        Actions:
        --------
        - Constructs a new filename using DATE-OBS, FILTER, OBJECT, and CCD ID.
        - Applies consistent formatting is:
        'Calib-{ccd}-{OBJECT}-{UTDATE}-{UTSTART}-{FILTER}-{EXPTIME}.fits'
        - Renames files on disk and updates any internal references.

        Output:
        -------
        - Files on disk are renamed to match the standard format.
        - Useful for organizing files and avoiding duplicates across nights/CCDs.

        Notes:
        ------
        - Should be used with caution if filenames are already used in another system.
        """

        if imlist_name == None:
            imlist = self.solved_images
            self.summary['CALNAME'] = self.summary['FILENAME'].copy()
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
                
                if self.summary['IMAGETYP'][idx].item() == 'Light Frame' :
                    self.summary['CALNAME'][idx] = newimage
                else:
                    self.summary['CALNAME'][idx] = inim
                
                print(f'Copy {inim} to {newimage}...')
                shutil.copy2(inim, newimage)
            print(f"Basic preprocessing of {self.ccd} is finished.")

            print('Remove previous check images...')
            os.system(f'rm {self.curdir}/Cal*.ch*.fits')
            
        elif imlist_name != None:
            imlist   = glob.glob(imlist_name)
            for inim in imlist :
                hdr     = fits.getheader(inim)
                EXPTIME = int(hdr['exptime'])
                UTDATE  = hdr['date-obs'][0:10]
                UTSTART = hdr['date-obs'][11:19]
                OBJECT  = hdr['OBJECT']
                FILTER  = hdr['FILTER']

                newimage = f'Calib-{self.ccd}-{OBJECT}-{str(UTDATE[0:4])}{str(UTDATE[5:7])}{str(UTDATE[8:10])}-{str(UTSTART[0:2])}{str(UTSTART[3:5])}{str(UTSTART[6:8])}-{FILTER}-{str(EXPTIME)}.fits'

                print(f'Copy {inim} to {newimage}...')
                shutil.copy2(inim, newimage)
            print(f"Basic preprocessing of {self.ccd} is finished.")  
        self.reduced_images = glob.glob('Cal*.fits') ; self.reduced_images.sort()

    def FluxScaling(self, imlist_name = None, zp0=25, zpkey='ZP'):
        """
        Normalize the flux of science images for consistent photometric scaling.

        Actions:
        --------
        - Computes the average flux level (e.g., background or source-based) for each image.
        - Scales images multiplicatively to match a common reference level.
        - Useful for aligning flux levels before stacking or image subtraction.

        Output:
        -------
        - Saves scaled images, possibly with a new suffix or overwrites.
        - Improves accuracy in subsequent steps like image subtraction or photometry.

        Notes:
        ------
        - Assumes relative photometry across images in the same filter and object.
        - Optional step depending on user's stacking or subtraction method.
        """
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
        Combine aligned science images into a single stacked image.

        Actions:
        --------
        - Groups images by filter and object.
        - Uses IRAF's `imcombine` or a similar method to median- or average-stack images.
        - Optionally applies rejection algorithms (e.g., `none`, `minmax`, etc.) based on self.reject.

        Output:
        -------
        - Produces a stacked FITS image for each group.
        - Improves signal-to-noise ratio and removes cosmic rays or artifacts.

        Notes:
        ------
        - Image alignment (WCS-based or pixel-based) should be performed beforehand.
        """

        if imlist_name == None:
            imlist      = self.fscaled_images
            objset  =  list(set(self.summary['OBJECT'][self.summary['IMAGETYP'] == 'Light Frame']))
            bandset =  list(set(self.summary['FILTER'][self.summary['IMAGETYP'] == 'Light Frame']))
            binset  =  list(set(self.summary['XBINNING'][self.summary['IMAGETYP'] == 'Light Frame']))
           

        elif imlist_name != None:
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
        Perform aperture photometry on science or stacked images.

        Actions:
        --------
        - Uses a configured photometry tool (e.g., SExtractor or custom script).
        - Extracts source fluxes, positions, and errors.
        - Optionally matches photometric results with known catalogs.

        Output:
        -------
        - Saves photometry catalogs for each input image.
        - Can be used for generating light curves or calibrating magnitudes.

        Notes:
        ------
        - Photometry configuration is loaded from `self.photpath`.
        - May depend on prior background subtraction or image quality checks.
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
                        phot(inim, target_catalog='', band=band, path=self.photpath, savecat=True, subprefix='hdCalib', ref='PS1', mykey='APER_1', sub=False, minarea=3, det_thresh=3, deb_nthresh=32, deb_mincont=0.01, backsize=128, backfiltersize = 5, backphoto_type='LOCAL',  ratio=1, onlyzp=True, snrcut=0.1, magup=12, maglow=17)

                    except:
                        try:
                            phot(inim, target_catalog='', band=band, path=self.photpath, savecat=True, subprefix='hdCalib', ref='APASS', mykey='APER_1', sub=False, minarea=3, det_thresh=3, deb_nthresh=32, deb_mincont=0.01, backsize=128, backfiltersize = 5, backphoto_type='LOCAL', ratio=1, onlyzp=True, snrcut=0.1, magup=11, maglow=15)
                        except:
                            print('No stars matched with APASS or PS1. PLEASE DO MANUALLY.')

                            # Ensure the target directory exists
                            os.makedirs('./photretry/', exist_ok=True)

                            # Get base name without .fits extension
                            basename = inim[:-5]

                            # Move files
                            shutil.move(inim, './photretry/')
                            shutil.move(f'{basename}.chaper.fits', './photretry/')
                            shutil.move(f'{basename}.chbkg.fits', './photretry/')
                            shutil.move(f'{basename}.chseg.fits', './photretry/')
                            shutil.move(f'{basename}.merge.cat', './photretry/')

                            print(f'{inim} is moved to photretry directory.')

                            self.reduced_images.remove(inim)
                else:
                    phot(inim, target_catalog='', band=band, path=self.photpath, savecat=True, subprefix='hdCalib', ref=ref, mykey='APER_1', sub=False, minarea=3, det_thresh=3, deb_nthresh=32, deb_mincont=0.01, backsize=128, backfiltersize = 5, backphoto_type='LOCAL',  ratio=1, onlyzp=True, snrcut=0.1, magup=magup, maglow=maglow)
                        
            else:
                print(f'{band}-band is not supported. TBD.')
        os.system('ls ./photretry/*.fits')
        print('Finished.')  

    def Subtract(self, imlist_name='Cal*.com.fits', template_dir=None):
        """
        Perform image subtraction using HOTPANTS.

        Actions:
        --------
        - Matches science images with template images from `self.template_dir`.
        - Applies PSF-matching and subtraction via HOTPANTS.
        - Produces residual images highlighting transient or variable sources.

        Output:
        -------
        - Saves subtracted FITS images to disk.
        - Optionally performs photometry on subtracted images.

        Raises:
        -------
        - RuntimeError if `self.template_dir` is None or does not contain valid templates.

        Notes:
        ------
        - Template images must be pre-aligned and photometrically matched.
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
        Archive processed images into a structured directory based on CCD and program type.

        Parameters:
        -----------
        imlist_name : str, optional
            If provided, only files matching this pattern will be archived.
            If not, the method uses self.reduced_images + self.combined_images + self.sub_images.

        Actions:
        --------
        - Reads image headers to determine object name, filter, and observing program.
        - Moves files to archive directories grouped by program (IMSNG, Transient, Exoplanet) or CCD type.
        - Sets appropriate permissions for access.

        Output:
        -------
        - Processed images are moved to `/mnt/dataset/obsdata/[Program or CCD]/[Object]/[CCD]/[Filter]/`.

        Notes:
        ------
        - Assumes the directory structure exists or is creatable.
        - Useful for centralizing data from multiple observing runs or instruments.
        """
        imsng_archive     = self.archive_paths['IMSNG']
        transient_archive = self.archive_paths['Transient']
        exoplanet_archive = self.archive_paths['Exoplanet']

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

            band    = hdr['FILTER']
            program = hdr.get('PROGRAM', '')

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
                if self.ccd.startswith('MAAO'):
                    archive_path = os.path.join(self.archive_paths['MAAO'], obj, self.ccd, band)
                elif self.ccd.startswith('PNUO'):
                    archive_path = os.path.join(self.archive_paths['PNUO'], obj, self.ccd, band)
                else:
                    raise ValueError(f"No archive path defined for CCD '{self.ccd}' and program '{program}'")
            os.makedirs(archive_path, exist_ok=True)
            dest_file_path = os.path.join(archive_path, os.path.basename(inim))

            shutil.copy2(inim, dest_file_path)

            os.system(f'chmod -R 777 {archive_path}')
        print('Finished.')

def main(ccd='PNUO_C361K'):
    import os
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
