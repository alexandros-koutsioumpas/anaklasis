# anaklasis 

[_anaklasis_](https://github.com/alexandros-koutsioumpas/anaklasis) is a set of open-source _Python3_ scripts (with _fortran90_ extensions) that facilitate a range of specular neutron and x-ray reflectivity calculations, involving the generation of theoretical curves and the comparison/fit of interfacial model reflectivity against experimental datasets.  The _ref_ module, contains three callable functions, _ref.calculate_ for generating theoretical reflectivity curves, _ref.compare_ for comparison of experimental data with theoretical curves and _ref.fit_ for refinement of experimental data against a defined model. Execution takes place by defining the interfacial model and instrumental parameters as lists in a simple _Python_ script and by passing them as arguments to the desired function.

In the examples folder many scripts with calculations and refinements performed by the program can be found. Addiitonaly some _Jupyter notebooks_ explaining the input and output are included.

Full description of used methods will be reported in the form of an open-access article.

## Installation

**Linux**

- Install _Python_ >= 3.7 and _gfortran_
- Install NumPy
- then download _anaklasis_ and install thought the terminal

```bash
python3 setup.py install --user
```

**macOS**

- Install _python_ >= 3.7 from [python.org](https://www.python.org/downloads/)
- Install _gfortran_ compiler. An easy way is to use the installers provided by _fxcoudert_ at [github](https://github.com/fxcoudert/gfortran-for-macOS)
- Install NumPy  

```bash
python3 -m pip install numpy
```

- then download _anaklasis_, navigate to the proper folder  and install thought the terminal

```bash
python3 setup.py install
```

Note that if you prefer you may use [MacPorts](https://www.macports.org) or [Homebrew](https://brew.sh) for the installation of _gfortran_ compiler.

**Windows 10 (without installing _gfortran_ )**

For convenience a _wheel_ with a pre-compiled _fortran_ extension is provided (folder `\win_wheel`) for _Python_ version 3.9 (more _Python_ versions to come). For installing it follow the steps below:

- Install _python_ 3.9 from [python.org](https://www.python.org/downloads/) and do not forget to include _python_ in the systen path.
- Install NumPy and upgrade setuptools
```bash
py -m pip install --upgrade pip setuptools

py -m pip install numpy
```

then navigate in the `\win_wheel` folder and install through the command prompt

```bash
py -m pip install anaklasis-1.5.0-cp39-cp39m-win_amd64.whl
```

In case you prefer *Anaconda* on *Windows* just make sure you have _setuptools_ and _NumPy_ on the system and just install the wheel.

**Windows (with _gfortran_ installation)**

- Install _python_ >= 3.7 from [python.org](https://www.python.org/downloads/) and do not forget to include _python_ in the systen path.
- Install NumPy and upgrade setuptools
```bash
py -m pip install --upgrade pip setuptools

py -m pip install numpy
```
- Install _MSYS2_ following these [instructions](https://www.scivision.dev/install-msys2-windows). Make sure that you include _c:\msys64\mingw64\bin_ in the system PATH (edit environement variables in control panel) and that you get the latest packages by using the appropriate commands (_pacman_ _-Syuu_) in the _MSYS2_ console.
- After _MSYS2_ is installed and all packages updated, install _gfortran_ from the _MSYS2_ console

```bash
pacman -S mingw64/mingw-w64-x86_64-gcc-fortran
```

- Finally download _anaklasis_, navigate to the proper folder and install thought the command prompt

```bash
py setup.py install build_ext --compiler=mingw32
```

In case you prefer *Anaconda* on *Windows* [here are instructions](https://python-at-risoe.pages.windenergy.dtu.dk/compiling-on-windows/configuration.html) for setting up the _fortran_ compiler. 

**Windows 10 (using WSL)**

An additional option for Windows 10 users is to use the Windows Subsystem for Linux (WSL), install a Linux distribution (like Ubuntu) and follow the installation instructions presented above for Linux.

## Running _anaklasis_ in the cloud

You may use the following **Binder** link in order to perform _anaklasis_ calculations inside _jupyter notebooks_ that run on the cloud .

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/alexandros-koutsioumpas/anaklasis/HEAD?filepath=templates_binder)

You may find templates for calculations and refinements in the form of _jupyter notebooks_ that you can modify according to your needs. You may also upload data files in order to use them in your refinements. Resources on **Binder** are a bit limited but you are able to experiment with _anaklasis_ without performing a local installation on your machine.

## Getting help

In the `/docs` folder you may find the program's API in various formats.

Also after installing *anaklasis*, you may install [pdoc](https://pdoc3.github.io/pdoc/) and from the console run

```python
pdoc anaklasis
```

that will open a html version of the API in your browser.

In this [YouTube link](https://www.youtube.com/watch?v=ieulImJUK5o) you may also find a video tutorial for fitting the [ORSO example data](https://github.com/reflectivity/reflectivity.github.io/blob/master/workshops/workshop_2021/ORSO_example.ort) with _anaklasis_.

## Templates

In the directory `/script_templates` you may find commented scripts that may guide you for writing _anaklasis_ analysis scripts.

## Post-installation tests

In the folder `/tests` you will find a _Python_ script `tests.py` that runs tests of the core calculations used by all three _anaklasis_ functions and reports if the results are the same as those obtained by version 1.3 of the package. You may run the tests using

```bash
python3 tests.py
```

If the run is succesful you will receive a `All anaklasis tests passed!` message at the end.

## Note about running Bootstrap analysis

On _macOS_ and _Linux_ before performing a bootstrap analysis you might need to increase the user resource limit using the command

```bash
ulimit -n 2048
```

## License
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)
