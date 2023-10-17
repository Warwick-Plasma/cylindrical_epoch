import os
import sys
from distutils.sysconfig import get_python_lib
from distutils.core import setup, Extension
import subprocess as sp
try:
    import numpy
    gotnumpy = True
except:
    gotnumpy = False


def get_numpy_dir():
    if gotnumpy:
        numpydir = os.path.dirname(numpy.__file__)
        incdir = os.path.join(numpydir, 'core', 'include')
        for r, d, fl in os.walk(incdir):
            for f in fl:
                if f == 'arrayobject.h':
                    return os.path.realpath(os.path.join(r, '..'))
        for r, d, fl in os.walk(numpydir):
            for f in fl:
                if f == 'arrayobject.h':
                    return os.path.realpath(os.path.join(r, '..'))
    sys.path.insert(0, get_python_lib(standard_lib=1))
    for path in sys.path:
        for r, d, fl in os.walk(path):
            for f in fl:
                if f == 'arrayobject.h':
                    return os.path.realpath(os.path.join(r, '..'))
    print('Unable to build python module. Numpy directory not found.')
    sys.exit(1)


sdfdir = os.path.join('..', 'C')

srcfiles = ['sdf_python.c']

incdirs = [get_numpy_dir()] + [os.path.join(sdfdir, 'src')]
libdirs = [os.path.join(sdfdir, 'lib')]


def get_version():
    commit_id = "UNKNOWN"
    version = "UNKNOWN"

    with open("commit_info.h") as f:
        lines = f.readlines()
        for l in lines:
            if l.find('SDF_COMMIT_ID') != -1:
                commit_id = l.split('"')[1]
                version = commit_id.lstrip('v').split('-')[0]
            elif l.find('SDF_COMMIT_DATE') != -1:
                commit_date = l.split('"')[1]

    with open("sdf_helper/_version.py", "w") as f:
        f.write('__version__ = "{}"\n'.format(version))
        f.write('__commit_id__ = "{}"\n'.format(commit_id))
        f.write('__commit_date__ = "{}"\n'.format(commit_date))

    return commit_id

get_version()

setup(name="sdf", version="1.0",
      ext_modules=[Extension("sdf", srcfiles, include_dirs=incdirs,
                   library_dirs=libdirs, libraries=['sdfc'])],
      packages=["sdf_helper"], py_modules=["sdf_legacy"])
