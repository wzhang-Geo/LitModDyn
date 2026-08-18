# LitModDyn

**_LitModDyn_** is a Python-based post-processing tool designed to compute present-day mantle flow and dynamic topography. It incorporates rheological layering derived from integrated geophysical–petrological modelling ([_LitMod2D_2.0_](https://github.com/ajay6763/LitMod2D_2.0_package_dist_users)), ensuring consistency between density, temperature, seismic structure, and viscosity. This allows for a more realistic representation of lithosphere–mantle coupling and its impact on surface topography.

The code solves the steady-state Stokes equations in a two-dimensional Cartesian domain using finite-difference and marker-in-cell techniques. It assumes an incompressible, viscous mantle and neglects inertial forces, which is appropriate for long-term mantle flow at low Reynolds numbers.

## Repository structure
The main files and directories are organized as follows:
```text
LitModDyn/
│
├── README.md
├── LICENSE
│
├── main.py
├── main_lit.py
├── Solver.py
├── Interpolation.py
├── Post_iso_dyn.py
│
├── Zhang_Northern_ext.py
├── Zhang_Southern_ext.py
│
├── input/
│   ├── Zhang_etal_2022_NorthPro/
│   │   ├── Best_Model/
│   │   └── Best_Model_No_Slabs/
│   │
│   └── Zhang_etal_2024_SouthPro/
│       ├── Best_Model/
│       └── Best_Model_No_Slabs/
│
├── Benchmark_Fallingblock/
│
├── Zhang_etal_2022_NorthPro/
└── Zhang_etal_2024_SouthPro/
```
---
## Main source files
`main.py` – main LitModDyn program used to calculate the present-day mantle-flow field.

`Solver.py` – finite-difference solver for the two-dimensional Stokes equations.

`Interpolation.py` – marker-in-cell interpolation routines between markers and the computational grid.

`main_lit.py` – rheological functions and auxiliary routines used in the mantle-flow calculation.

`Zhang_Northern_ext.py` – initialization and model configuration for the Northern profile.

`Zhang_Southern_ext.py` – initialization and model configuration for the Southern profile.

`Post_iso_dyn.py` – post-processing routines used to calculate dynamic topography from the LitModDyn results.


## Input files
The `input/` directory contains the model files required to initialize the principal simulations presented in the manuscript.
```text
input/
├── Zhang_etal_2022_NorthPro/
│   ├── Best_Model/
│   └── Best_Model_No_Slabs/
│
└── Zhang_etal_2024_SouthPro/
    ├── Best_Model/
    └── Best_Model_No_Slabs/
```
The Northern profile is based on the geophysical–petrological model of  [_Zhang et al. (2022)_](https://doi.org/10.1029/2022JB024800), whereas the Southern profile is based on  ([_Zhang et al. (2024)_](https://doi.org/10.1029/2023JB028435)).


For each profile:
- `Best_Model/` contains the input files for the best fitting model.
- `Best_Model_No_Slabs/` contains the corresponding model without mantle anomalies. 


## Installation
LitModDyn can be downloaded from the archived Zenodo release or cloned from the GitHub repository.
The code is compatible with Python 3 and requires the following Python libraries:
```text
numpy
scipy
matplotlib
multiprocessing
```
## Reproducibility workflow
To reproduce the principal results, an independent user should therefore:
- Clone or download the archived LitModDyn release.
- Create a Python 3 environment and install the required external libraries.
- Confirm that the complete input datasets are present in `input/`.


## Citing
If you use this tool please cite:

Zhang, W., Jiménez-Munt, I., Negredo, A. M., García-Castellanos, D., Ortega Gelabert, O., Vergés, J., Sharma, M., & Torné, M. (2026).  
LitModDyn: Quantifying dynamic topography induced by slab–mantle interaction beneath the Adria microplate. *Frontiers in Earth Science*. **(in revision)**

Zhang, W., Jiménez-Munt, I., Negredo, A. M., García-Castellanos, D., Ortega Gelabert, O., Vergés, J., Sharma, M., & Torné, M. (2026).  
*LitModDyn: Quantifying dynamic topography induced by slab–mantle interaction*. Zenodo. https://doi.org/10.5281/zenodo.20746376

