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
    eff1 = []
    u_eff1 = []
    eff2 = []
    u_eff2 = []
    n_ph = 1E6
    z_range = np.arange(2,76,2) # z positions from 0 to 60cm in steps of 2.5cm
    for z in z_range:
        generator = DefaultPhotonGenerator(n_ph, 0, 0, 6, z)
        generator.generate()

        #WOM geometry values in cm or deg, according to Upgrade geometry
        #length, a_inner, b_inner, alpha_inner, x_inner, y_inner, a_outer, b_outer, alpha_outer
        model = EllipticModel(76, 5.5, 5.5, 0, 0, 0, 5.75, 5.75, 0)
        #n1 refr index WOM, n2 refr index environment (air)
        #model, n1, n2, lambda_abs, lambda_sc, generator
        simulation = Simulation(model, 1.46, 1.33, 416, 135, generator)
        result = simulation.simulateData()
        # number of cathode hits = sum(status)
	# status = 1 or result[2]==1, upper end of tube
        n_hit1 = np.where(result[2]==1,1,0).sum()
	# status = 5 or result[2]==5, lower end of tube
        n_hit2 = np.where(result[2]==5,1,0).sum()
	# efficiency = cathod hits / all photons
        eff1.append(n_hit1/n_ph)
        eff2.append(n_hit2/n_ph)
	# uncertainty on efficiency, total number of photons fix
        u_eff1.append(np.sqrt(n_hit1)/n_ph)
        u_eff2.append(np.sqrt(n_hit2)/n_ph)
    eff1 = np.array(eff1)
    u_eff1 = np.array(u_eff1)
    eff2 = np.array(eff2)
    u_eff2 = np.array(u_eff2)
    file = open("./files/sim_eff.pkl", "wb")
    data = [z_range, eff1, u_eff1, eff2, u_eff2]
    pickle.dump(data, file)
    file.close()


if __name__ == '__main__':
    main()
