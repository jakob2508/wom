# MC simulation
Underlying python-based MC written by Anna Pollmann and tested by Yuriy Popovych. Repository can be found on GitLab if access is provided
(https://gitlab.rlp.net/WOM/transmission-ray-tracing). The simulation is used to simulate the angular response of the sensor by shining light on the sensor from
a plane with different zenith angles. The simulation takes transmission losses and Fresnel effects into account. The simulation also takes coupling efficiencies
into account which are not considered in this setup because we later renormalize the obtained angular sensitivity curve.

# sim.py
Has to be placed in the same folder where Optics.py, Material.py, etc. are defined. With 1 deg zenith steps and 1E5 photons in each step the simulation takes
about 30 min. Saves the zenith angle, efficiency and uncertainty on the efficiency in a dictionary in the file efficiency_vs_zenith.csv.

# plotter.py
Plots the angular sensitivity in theta and cos(theta) and fits a polinomial through the data. The coefficients of the fit are saved in as_140.0.

# To Do
- include sampling from custom angular sensitivity file in ppc
- so far sin(theta) is assumed which is a good approximation but Fresnel effects quick in for small angles
