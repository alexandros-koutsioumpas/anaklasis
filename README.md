# anaklasis 

_anaklasis_ is a set of open-source _python3_ scripts (with _fortran90_ extensions) that facilitate a range of specular neutron and x-ray reflectivity calculations, involving the generation of theoretical curves and the comparison/fit of interfacial model reflectivity against experimental datasets.  The _anaklasis_ module, contains three callable functions, _anaklasis.calculate_ for generating theoretical reflectivity curves, _anaklasis.compare_ for comparison of experimental data with theoretical curves and _anaklasis.fit_ for refinement of experimental data against a defined model. Execution takes place by defining the interfacial model and instrumental parameters as lists in a simple _python_ script and by passing them as arguments to the desired function.

In the examples folder many scripts with calculations and refinements performed by the program can be found. Addiitonaly some _Jupyter notebooks_ explaining the input and output are included.

Full description of used methods will be reported in the form of an open-access article.

## Installation

**Linux**

- Install _Python_ > 3.7 and _gfortran_
- Install NumPy
- then download _anaklasis_ and install thought the terminal

```bash
python3 setup.py install --user
```

**macOS**

- Install _python_ > 3.7 from [python.org](https://www.python.org/downloads/)
- Install _gfortran_ compiler. An easy way is to use the installers provided by **fxcoudert** at [github](https://github.com/fxcoudert/gfortran-for-macOS)
- Install NumPy  

```bash
python3 -m pip install numpy
```

- then download _anaklasis_, navigate to the proper folder  and install thought the terminal

```bash
python3 setup.py install
```

**Windows**

- Install _python_ > 3.7 from [python.org](https://www.python.org/downloads/) and do not forget to include _python_ in the systen path.
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


## Getting help

After installing the program, you may open a Python shell and import anaklasis

```python
from anaklasis import anaklasis
```

and then use the help command for the function you are interested in, like for example:

```python
help(anaklasis.fit)
```

A description of required input and output will be displayed. A file with the program's API is also provided (_anaklasis-API.txt_).

## Note about running Bootstrap analysis

On _macOS_ and _Linux_ before performing a bootstrap analysis you might need to increase the user resource limit using the command

```bash
ulimit -n 2048
```

## License
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)
