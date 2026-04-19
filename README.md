# LitModDyn

**_LitModDyn_** is a Python-based post-processing tool designed to compute present-day mantle flow and dynamic topography. **LitModDyn** incorporates rheological layering derived from integrated geophysical–petrological modelling ([_LitMod2D_2.0_](https://github.com/ajay6763/LitMod2D_2.0_package_dist_users)), ensuring consistency between density, temperature, seismic structure, and viscosity. This allows for a more realistic representation of lithosphere–mantle coupling and its impact on surface topography.

The code solves the steady-state Stokes equations in a two-dimensional Cartesian domain using finite-difference and marker-in-cell techniques. It assumes an incompressible, viscous mantle and neglects inertial forces, which is appropriate for long-term mantle flow at low Reynolds numbers.

