from setuptools import setup, find_packages

setup(
    name='lgpyred',
    version='0.1',
    description='LGPY REDuction pipeline',
    author='Gu Lim',
    packages=find_packages(),
    include_package_data=True,
	install_requires=[
		'numpy==1.26.4',
		'astropy==6.0.1',
		'pyraf==2.2.2',
		'paramiko==3.4.0',
		'matplotlib==3.8.4',
		'scipy==1.12.0',
		'astroquery==0.4.7',
		'pillow==9.2.0',
	],
    entry_points={
        'console_scripts': [
            'lgpyred = lgpyred.lgpyred:main'
        ]
    },
)

