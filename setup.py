
import setuptools 
from numpy.distutils.core import setup, Extension


ext = [Extension(name='fortran_ref',
                 sources=['anaklasis/fortran_ref_functions.f90'],
                 f2py_options=['--quiet'],
                 extra_f90_compile_args=['-w'])]

setup(
	name='anaklasis',
	version='1.4.6',
	author='Alexandros Koutsioumpas',
	author_email='a.koutsioumpas@fz-juelich.de',
	description='Neutron and X-ray reflectivity calculations',
	packages=setuptools.find_packages(),
	ext_modules=ext,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
       ],
    python_requires='>=3.7',
    install_requires=['numpy>=1.18','scipy>=1.4','matplotlib','numdifftools>=0.9.39','sympy>=1.6.2','emcee>=3.0','tqdm','corner'],
	)



