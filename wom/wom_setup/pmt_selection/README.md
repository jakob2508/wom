# ppc_wom/wom_setup/pmt_selection

# ALGO
This repository contains the simulation script sim.py which runs the ALGO simulation returning the position, direction and time of photons traversing the inner (wavelength-shifting) tube of the WOM.
ALGO neither simulates how the photons enter the tube nor the wavelength-shift but assumes instead an isotropic emission point inside the tube. ALGO is optimized for GPU use and reaches a throughput of O(1E7) photons per second.
Thus ALGO only consideres relfection absorption and scattering inside the tube. As transmission losses, PMT effciencies are already covered by the effective area parametrization we can study the effects of PMT selection
probability as a function of hit poisition on the inner cylinder as well as the time distribution.

Further information on ALGO can be found in this master thesis by Florian Thomas: https://drive.google.com/drive/u/0/folders/0Bww7VbdYb87BWWloZmpCd0pJME0?resourcekey=0-MtNXpLTW2ntITKZUrD8P8A
In addition if access is given there is also GitLab repository (https://gitlab.rlp.net/WOM/algo) for further read.

The simulation was run on a GPU on the Mainz computing cluster. For information on how to set up ALGO on a GPU please contact the WOM group.

# sim.py
The following geometry has been used: r_inner = 5.5 cm, r_outer = 5.75 cm, length = 76 cm, n_quartz = 1.46, n_filling = 1.33
Scattering length: l_sca = 135 cm
Absorption length: l_abs = 416 cm
As given by Yuriy's latest combined fit of the attenuation length.
Loops over different cylinder heights z on the inner tube and saves the efficiency (number of hits/number of sim. photons) for both tube ends and it's uncertainty.

# plotter.py
Plots the efficiency as a function of the cylinder position z. Renormalizes the efficiency such that the total is 100% and fits a linear function through the data. The parameters of that function are saved in
ps_140.0 and are used in ppc to select a PMT given the z hit position on the cylinder by using rejection sampling. 

# files/sim_eff.pkl
Contains the results of sim.py:

z_position eff_end1 u_eff_end1 eff_end2 u_eff_end2

