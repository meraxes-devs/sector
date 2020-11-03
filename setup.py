from distutils.core import setup
from setuptools import find_packages, Extension
from Cython.Distutils import build_ext


exts = [Extension("sector.utils",
                  sources = ["sector/utils.pyx"],
                  include_dirs = ["sector/"]),
        Extension("sector.sfh",
                  sources = ["sector/sfh.pyx"],
                  include_dirs = ["sector/"]),
        Extension("sector.sector",
                  sources = ["sector/clib/tools.c",
                             "sector/clib/sector.c",
                             "sector/clib/spectra.c",
                             "sector/clib/dust.c",
                             "sector/clib/igm.c",
                             "sector/sector.pyx"],
                  include_dirs = ["sector/"],
                  extra_compile_args = ['-fopenmp', '-lhdf5', '-lhdf5_hl'],
                  extra_link_args = ['-fopenmp', '-lhdf5', '-lhdf5_hl']),
        Extension("sector.dust.dust",
                  sources = ["sector/clib/tools.c",
                             "sector/clib/sector.c",
                             "sector/clib/spectra.c",
                             "sector/clib/dust.c",
                             "sector/clib/igm.c",
                             "sector/dust/dust.pyx"],
                  include_dirs = ["sector/"],
                  extra_compile_args = ['-fopenmp', '-lhdf5', '-lhdf5_hl'],
                  extra_link_args = ['-fopenmp', '-lhdf5', '-lhdf5_hl'])]

for e in exts:
    e.cython_directives = {"embedsignature":True, "language_level": "2"}

exec(open('sector/version.py', 'r').read())

setup(
    name = 'sector',
    version = __version__,
    description = 'Package to compute SEDs for semi-analytic models',
    author = 'Yisheng Qiu',
    author_email = 'yishengq@student.unimelb.edu.au',
    ext_modules = exts,
    cmdclass = {'build_ext': build_ext},
    packages = find_packages(),
    package_data = {'sector/filters':['*.npy'], 'sector/sedlib/':['*.hdf5']}
)
