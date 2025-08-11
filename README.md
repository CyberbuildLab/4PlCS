# 4PlCS
This repository contains implementations of the method presented in "4-Plane Congruent Sets for Automatic Registration of As-Is 3D Point Clouds with 3D BIM Models" (https://doi.org/10.1016/j.autcon.2018.01.014) for registrating a dense point cloud (e.g. from last scanning or photogrammetry) with a 3D BIM model of 
More specifically, it includes:
- The original Matlab code used to produce the results in the manuscript. This code is provided "raw" and was not commented well. This assumes that the input 3D BIM model is converted in OBJ format.
- A Jupyter Notebook implementing that method too, but assuming that the input 3D BIM model is in IFC format. THe Jupyter Notebook is annotated to explain each step and the parameters involved. IMPORTANTLY: the Jupyter Notebook implementation was developed independently from the Matlab implementation; both may thus differ in their implementation (e.g. assumptions about the pcd normal vectors)

Any feedback, suggestions or even contributions (e.g. Python library implementation, or structuring of the Matlab code) would be welcome. Don't hesitate to contact me.

## Acknowledgement
```
@article{BUENO2018120,
author = {Martín Bueno and Frédéric Bosché and Higinio González-Jorge and Joaquín Martínez-Sánchez and Pedro Arias},
title = {4-Plane congruent sets for automatic registration of as-is 3D point clouds with 3D BIM models},
journal = {Automation in Construction},
volume = {89},
pages = {120-134},
year = {2018},
issn = {0926-5805},
doi = {https://doi.org/10.1016/j.autcon.2018.01.014}}
```
