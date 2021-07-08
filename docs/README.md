# anaklasis 

[_anaklasis_](https://github.com/alexandros-koutsioumpas/anaklasis) is a set of open-source _Python3_ scripts (with _fortran90_ extensions) that facilitate a range of specular neutron and x-ray reflectivity calculations, involving the generation of theoretical curves and the comparison/fit of interfacial model reflectivity against experimental datasets.  The _ref_ module, contains three callable functions, _ref.calculate_ for generating theoretical reflectivity curves, _ref.compare_ for comparison of experimental data with theoretical curves and _ref.fit_ for refinement of experimental data against a defined model. Execution takes place by defining the interfacial model and instrumental parameters as lists in a simple _Python_ script and by passing them as arguments to the desired function.


In the [project repository page](https://github.com/alexandros-koutsioumpas/anaklasis) you may find instructions for installing the package on all major platforms (_Linux_, _macOS_, _Windows_) and also documentation and examples in the form of scripts and _Jupyter_ notebooks.

As an alternative you may use the following **Binder** link in order to perform _anaklasis_ calculations inside _jupyter notebooks_ that run on the cloud .

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/alexandros-koutsioumpas/anaklasis/HEAD?filepath=templates_binder)

Full description of used methods will be reported in the form of an open-access article.

In this [YouTube link](https://www.youtube.com/watch?v=ieulImJUK5o) you may also find a video tutorial for fitting the [ORSO example data](https://github.com/reflectivity/reflectivity.github.io/blob/master/workshops/workshop_2021/ORSO_example.ort) with _anaklasis_.


