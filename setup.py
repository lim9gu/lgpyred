from setuptools import setup, find_packages

setup(
    name='lgpyred',
    version='0.1',
    description='LGPY Image Reduction Pipeline',
    author='Gu Lim',
    packages=find_packages(),
    include_package_data=True,
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

