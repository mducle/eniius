#!/usr/bin/env

from setuptools import setup

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

config = dict(
    name='eniius',
    version=versioneer.get_version(),
    author='Duc Le',
    author_email='duc.le@stfc.ac.uk',
    description='A utility for embedding neutron instrument information using (nx)spe files.',
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    packages=['eniius'],
    install_requires = ['nexusformat>=0.7.8', 'mcstasscript>=0.0.54'],
    extras_require = {},
    cmdclass=cmdclass,
    entry_points={'console_scripts': ['eniius = eniius:main']},
    url="https://github.com/mducle/eniius",
    zip_safe=False,
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Microsoft :: Windows :: Windows 10",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics",
    ]
}

try:
    setup(**config)
except CalledProcessError:
    print("Failed to build the extension!")
