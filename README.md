# lgpyred

**lgpyred** is a pipeline for reducing astronomical imaging data, especially for PNUO, MAAO, etc. It is designed to work with a combination of powerful external tools such as **IRAF**, **astrometry.net**, **wcsremap**, and **hotpants (optional)**, in addition to Python-based libraries like `astropy` and `numpy`.  

This pipeline was developed to automate standard procedures in observational astronomy, including astrometric calibration, image alignment, subtraction, and photometry.

<br>

Developed by **Gu Lim**  
Dept. of Earth Sciences, Pusan National University


## Installation Guide

To use **lgpyred**, follow the steps below to install the required Python packages and prepare your environment.

### 1. Clone the repository and install the Python package

If you have not already done so, clone the repository using Git and install the package locally using `pip`.

```bash
git clone https://github.com/lim9gu/lgpyred.git
cd lgpyred
pip install .
```

This will install `lgpyred` along with the following required Python dependencies:

- `astropy==6.0.1` – Core astronomy library for Python  
- `astroquery==0.4.7` – Tools to access online astronomical databases  
- `matplotlib==3.8.4` – Plotting library for visualization  
- `numpy==1.26.4` – Fundamental package for numerical computations  
- `scipy==1.12.0` – Scientific computing tools  
- `pyraf==2.1.13` – IRAF functionality in Python (deprecated but still used)  
- `paramiko==3.4.0` – SSH2 protocol library for remote operations  
- `pillow==9.2.0` – Python Imaging Library fork, for handling images

> **Note:** It is strongly recommended to use a Python virtual environment (e.g., via `venv` or `conda`) before installation to avoid conflicts with other Python packages.

### 2. Install required external software

The pipeline also depends on several standalone external tools that are not included in the Python package. These must be installed separately on your system:

- `IRAF 2.17 (https://github.com/iraf-community)
- `Astrometry.net (Lang et al. 2010)
- `wcsremap (by Andrew Becker)
- `hotpants (by Andrew Becker)

Please refer to the docs/ folder for detailed step-by-step installation guides for each tool.

### 3. Installation

Clone this repository and install using pip:

```bash
git clone https://github.com/lim9gu/lgpyred.git
cd lgpyred
pip install .
```

This will install lgpyred along with its required dependencies...

### 4. Example usage
Once installed, you can use lgpyred as a command-line tool:

```python
lgpyred --input image.fits --mode reduce --date 20240625
```
Alternatively, you can also use it as a Python module:

```python
from lgpyred.lgpyred import Red

red = Red(ccd='PNUO_C361K')
main()
```

### 5. Supported Platforms
This pipeline has been tested on Ubuntu Linux 24 LTS.
Support for macOS or Windows may require additional adjustments, especially for installing external tools like IRAF or Hotpants.

License
Distributed under the MIT License. See LICENSE file for details.
