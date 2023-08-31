
import setuptools 
from numpy.distutils.core import setup, Extension

try:
    ext = [Extension(name='fortran_ref',
                     sources=['anaklasis/fortran_ref_functions.f90'],
                     f2py_options=['--quiet'],
                     extra_f90_compile_args=['-w'])]

    setup(
    	name='anaklasis',
    	version='1.6.0',
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
        zip_safe=False,
        install_requires=['numpy<=1.22.4','scipy>=1.4','matplotlib','numdifftools>=0.9.39','sympy>=1.6.2','emcee>=3.0','tqdm','corner'],
    	)

    print(' ')
    print('anaklasis was installed and the FORTRAN extensions for reflectivity calculations were successfully built!')
except:

    setup(
        name='anaklasis',
        version='1.6.0',
        author='Alexandros Koutsioumpas',
        author_email='a.koutsioumpas@fz-juelich.de',
        description='Neutron and X-ray reflectivity calculations',
        packages=setuptools.find_packages(),
        classifiers=[
            "Programming Language :: Python :: 3",
            "Operating System :: OS Independent",
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
           ],
        python_requires='>=3.7',
        zip_safe=False,
        install_requires=['numpy<=1.22.4','scipy>=1.4','matplotlib','numdifftools>=0.9.39','sympy>=1.6.2','emcee>=3.0','tqdm','corner'],
        )

    print('*** Warning ***')
    print(' ')
    print('no FORTRAN compiler was detected and anaklasis was installed without building the FORTRAN extensions for reflectivity calculations..')
    print(' Install the Numba package to accelerate the calculation engine, otherwise the package will still run but with a very slow pure Python calculation engine.')


