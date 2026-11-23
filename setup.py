
import os
import shutil
import subprocess
import sys
import sysconfig

import setuptools
from setuptools import setup
from setuptools.command.build_py import build_py as _build_py
from setuptools.dist import Distribution

HERE = os.path.abspath(os.path.dirname(__file__))
# Must be absolute: f2py's meson backend copies sources into its own
# scratch build directory, and a relative path does not survive that.
F90_SOURCE = os.path.join(HERE, 'anaklasis', 'fortran_ref_functions.f90')
MODULE_NAME = 'fortran_ref'

# anaklasis/ref.py does `from fortran_ref import ...` (an absolute, top-level
# import), so the compiled extension has to end up as a standalone top-level
# module next to the anaklasis package, not inside it.


def _compiled_module_name():
    return MODULE_NAME + (sysconfig.get_config_var('EXT_SUFFIX') or '.so')


def build_fortran_extension():
    """Compile the FORTRAN reflectivity engine with f2py.

    numpy.distutils (and distutils itself) is gone, so the extension is
    built by shelling out to f2py's own meson/ninja backend instead of
    going through numpy.distutils.core.Extension.
    """
    try:
        subprocess.run(
            [sys.executable, '-m', 'numpy.f2py', '-c', F90_SOURCE,
             '-m', MODULE_NAME, '--quiet'],
            cwd=HERE, check=True,
        )
    except Exception:
        return False
    return os.path.exists(os.path.join(HERE, _compiled_module_name()))


class BinaryDistribution(Distribution):
    """Tells setuptools/wheel this package contains a compiled extension
    (built out-of-band by f2py), so wheels get a platform-specific tag
    instead of the default pure-Python 'py3-none-any'."""

    def has_ext_modules(self):
        return True


class build_py(_build_py):
    def run(self):
        built = build_fortran_extension()
        super().run()

        if built:
            shutil.copy2(
                os.path.join(HERE, _compiled_module_name()),
                self.build_lib,
            )
            print(' ')
            print('anaklasis was installed and the FORTRAN extensions for reflectivity calculations were successfully built!')
        else:
            print('*** Warning ***')
            print(' ')
            print('no FORTRAN compiler was detected and anaklasis was installed without building the FORTRAN extensions for reflectivity calculations..')
            print(' Install the Numba package to accelerate the calculation engine, otherwise the package will still run but with a very slow pure Python calculation engine.')


setup(
    name='anaklasis',
    version='1.6.0',
    author='Alexandros Koutsioumpas',
    author_email='a.koutsioumpas@fz-juelich.de',
    description='Neutron and X-ray reflectivity calculations',
    license='GPL-3.0-only',
    packages=setuptools.find_packages(),
    package_data={'anaklasis': ['fortran_ref_functions.f90']},
    cmdclass={'build_py': build_py},
    distclass=BinaryDistribution,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        'Environment :: Console',
        'Intended Audience :: Science/Research',
       ],
    python_requires='>=3.7',
    zip_safe=False,
    install_requires=['numpy>=1.22','scipy>=1.4','matplotlib','numdifftools>=0.9.39','sympy>=1.6.2','emcee>=3.0','tqdm','corner'],
    )
