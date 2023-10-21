#simulation
All files for the simulation can be found on the cobalt node under /home/jbeise/ppc/resources/ice/checks/eff_area.
The log file containing all hits is located in /home/jbeise/ppc/resources/ice/checks/eff_area/logs.

# ppc
This version of ppc has been modified such that absorption of the LED photons by the harness (ofla=-2), lower hemisphere (ofla=-3) and cables (ofla=-4) is disabled (see pro.cxx line 473-497). This was done in 
order to get the same results for negative and postive angles relative to the horizontal plane when measuring the angular sensitivity.
