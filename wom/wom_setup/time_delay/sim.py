#!/usr/bin/env python3

"""
Adaption of Ronja Schnurs example code that probes different z positions
on a prototype WOM tube (geomtry from WOM paper).
Data used to determine the z dependent detection efficiency as well as the
time delay for different z positions.

Run from ALGO root directory.
"""

__author__ = 'Jakob Beise'
__email__ = 'beisejak@physik.uni-mainz.de'

from pyalgo import (
    QuadraticSpline,
    QuadraticSplineModel,
    Simulation,
    PhotonLoader,
    DefaultPhotonGenerator,
    Intervall,
    EllipticModel
)
import math
import numpy as np
import time
import pickle
from tools import falke, parabola

def main():
    #use default photon generator - currently implements a point light source
    #implement a modified version in CUDA code if you want to modify the light source
    #Parameters are
    #nPhotons, sigma, lx, ly, lz, seed
    #{lx, ly, lz} defines light position, seed is optional
    #photon hit position z in cm, half of paint layer of thickness 17µm = 0.5µm
    times1 = []
    times2 = []
    n_ph = 1E6
    z_range = np.arange(2,76,2) # z positions from 0 to 60cm in steps of 2.5cm
    for z in z_range:
        generator = DefaultPhotonGenerator(n_ph, 0, 0, 6, z)
        generator.generate()

        #WOM geometry values in cm or deg
        #length, a_inner, b_inner, alpha_inner, x_inner, y_inner, a_outer, b_outer, alpha_outer
        model = EllipticModel(76, 5.5, 5.5, 0, 0, 0, 5.75, 5.75, 0)

        #n1 refr index WOM, n2 refr index environment (air)
        #model, n1, n2, lambda_abs, lambda_sc, generator
        simulation = Simulation(model, 1.46, 1.33, 416, 135, generator)
        result = simulation.simulateData()
	# number of cathode hits = sum(status)
        # status = 1 or result[2]==1, upper end of tube
        mask1 = result[2]==1
	# timing of photons that hit upper tube end
        time1 = result[1].T[2][mask1]
	# status = 5 or result[2]==5, lower end of tube
        mask2 = result[2]==5
	# timing of photons that hit lower tube end
        time2 = result[1].T[2][mask2]
        times1.append(time1)
        times2.append(time2)
    #times = np.array(times)
    #do stuff with result now
    #np.save('./files/sim_eff.npy', efficiency)
    file = open("./files/sim_delay.pkl", "wb")
    data = [z_range, times1, times2]
    pickle.dump(data, file)
    file.close()
    #result is tuple with
    #([x,y,z] <- np.array, [phi, theta, dt] <- np.array, exit_code <- np.array, detected <- int) ?
    #numpy.save('./files/sim_z={:.0f}cm_py'.format(z), numpy.array(result))


if __name__ == '__main__':
    main()
