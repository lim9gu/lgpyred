from setuptools import setup, find_packages

setup(
    name='lgpyred',
    version='0.1',
    description='LGPY Image Reduction Pipeline',
    author='Gu Lim',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'astropy',
        'pyraf',
        'paramiko',
    ],
    entry_points={
        'console_scripts': [
            'lgpyred = lgpyred.lgpyred:main'
        ]
    },
    include_package_data=True,
    package_data={
        'lgpyred.lgpytars': [
            'data/**/*.*',
            'astrom_config/*.*',
            'photconf/*.*',
            'sex_config/*.*',
            'reduction/*.*',
        ],
    },
    install_requires=[
        'numpy',
        'astropy',
        'pyraf',
        'paramiko',
    ],
    entry_points={
        'console_scripts': [
            'lgpyred = lgpyred.lgpyred:main'
        ]
    },
)
