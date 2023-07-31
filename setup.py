#!/usr/bin/env

import sys, os
import subprocess
import shutil
import versioneer
from setuptools import setup

def get_git_deps():
    # git clone dependent repos (McStas components and PyChop)
    cwd = os.getcwd()
    if not os.path.exists('McCode') or not os.path.exists(os.path.join('McCode', '.git')):
        rv = subprocess.run(['git', 'clone', 'https://github.com/McStasMcXtrace/McCode', '--depth=1'])
        if rv.returncode != 0:
            raise Exception(f'Could not clone McStas from github')
    comps_dst = os.path.join(cwd, 'eniius', 'mcstas-comps')
    if not os.path.exists(comps_dst):
        os.chdir('McCode')
        rv = subprocess.run(['git', 'restore', '--source=HEAD', ':/'])
        def copy_comp(src, dst, *, follow_symlinks=True):
            if src.endswith('.comp') and 'obsolete' not in src:
                shutil.copy2(src, dst, follow_symlinks=follow_symlinks)
        shutil.copytree('mcstas-comps', comps_dst, copy_function=copy_comp)
        for emptydir in ['data', 'examples', 'obsolete', 'share']:
            shutil.rmtree(os.path.join(comps_dst, emptydir))
        shutil.copytree(os.path.join('mcstas-comps', 'share'), os.path.join(comps_dst, 'share'))
        contrib_dst = os.path.join(comps_dst, 'contrib', 'ISIS_tables')
        for face in ['TS1_S01_Maps.mcstas', 'TS1_S04_Merlin.mcstas', 'TS2.imat']:
            shutil.copy2(os.path.join('mcstas-comps', 'contrib', 'ISIS_tables', face), contrib_dst)
        os.chdir(cwd)
    pychop_dst = os.path.join(cwd, 'eniius', 'pychop')
    if not os.path.exists(pychop_dst):
        if not os.path.exists('pychop-git'):
            rv = subprocess.run(['git', 'clone', 'https://github.com/mducle/pychop', 'pychop-git', '--depth=1'])
            if rv.returncode != 0:
                raise Exception(f'Could not clone PyChop from github')
        shutil.copytree(os.path.join('pychop-git', 'PyChop'), pychop_dst)


def recurse_data_files(rootdir, extn=None):
    datfiles = []
    for root, _, files in os.walk(rootdir):
        for ff in files:
            if extn is None or extn in ff:
                datfiles.append(os.path.join('..', root, ff))
    return datfiles


with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

config = dict(
    name='eniius',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author='Duc Le',
    author_email='duc.le@stfc.ac.uk',
    description='A utility for embedding neutron instrument information using (nx)spe files.',
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    packages=['eniius', 'eniius.pychop'],
    install_requires = ['nexusformat>=0.7.8', 'mcstasscript>=0.0.54', 'regex'],
    extras_require = {},
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
)

try:
    get_git_deps()
    data_files = recurse_data_files(os.path.join(os.path.dirname(__file__), 'eniius', 'mcstas-comps'))
    data_files += recurse_data_files(os.path.join(os.path.dirname(__file__), 'eniius', 'instruments'))
    data_files += recurse_data_files(os.path.join(os.path.dirname(__file__), 'eniius', 'pychop'), '.yaml')
    config['package_data'] = {'eniius': data_files}
    setup(**config)
except CalledProcessError:
    print("Failed to build the extension!")
