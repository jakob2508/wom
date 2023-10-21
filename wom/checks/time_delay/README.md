#simulation
To cross-check the implementation of the time delay we can use the same simulation results as for the angular sensitivity file ppc_wom/checks/ang_sens/sim/logs as that simulation covered the whole sensor. 
For the purpose of the time delay simulation the case theta = 90 deg would constitute a symmetric, somewhat uniform time distribution over the tube z position, where as theta = 10 deg should favour the upper 
PMT (PMT = 1) and theta = 170 deg the lower PMT (PMT = 0).

#analysis ./ana
Contains all scripts required to apply the time delay on the 'raw' hit files from PPC. WOM_time_delay.py contains all functions needed by check.py to perform some simple sanity checks. The plots are stored in 
/ana/plots
