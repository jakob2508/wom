import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import pandas as pd

from Optics import *

from Material import *
from Generator import *
from Object import *
from Propagator import *
from Cylinder import *

vessel_outer_dia=145+2*14
vessel_inner_dia=145
vessel_height=1300

tube_outer_dia=115
tube_inner_dia=115-2*2.5
tube_height=760

##########

n_air=1.0 # @ 589.29nm ?
n_quartz=1.46
n_coat=n_quartz # ?? check
n_glue=1.48
n_pmtglass=1.49 # ?? check with yurij
n_pfpe=1.32 # @ 589.29nm 
n_ice=1.33 # @ 589.29nm 

##########

number=1e5
plane_x=200
plane_y=1500
plane_dist=2000

###########

wls_eff=0.985002681771287
abs_scat=0.75
pmt_radial_acc=0.85
glue_transmission=0.9987637413055658
pmt_mean=0.20681970351875
all_factors=wls_eff*abs_scat*pmt_radial_acc*glue_transmission*pmt_mean

ice=Material("Ice", n_ice)
air=Material("Air", n_air)
glass=Material("Glass", n_quartz)
fill=Material("Fill", n_pfpe)# only temporarily, will be adjusted in loop

######################
### SETUP DETECTOR ###
######################

# Configure generator with rectangular plane
gen=Generator(number, "Flat", ("Rect", plane_x, plane_y,plane_dist))
environment=Environment("Ice",ice)
vessel=Cylinder("Vessel", environment, Vector(0,0,-0.5*vessel_height),glass,
                vessel_outer_dia/2 ,vessel_inner_dia/2, vessel_height,
               )

filling=Cylinder("Filling", vessel, Vector(0,0,80-0.5*vessel_height), fill,
                vessel_inner_dia/2, tube_outer_dia/2, tube_height,
                )
wls=Cylinder("WLS", filling, Vector(0,0,80-0.5*vessel_height), glass,
                tube_outer_dia/2,tube_inner_dia/2, tube_height,
                 "WLS"
            )

middle=Cylinder("Middle", wls, Vector(0,0,80-0.5*vessel_height), air,
                tube_inner_dia/2,0, tube_height, 
                )

listOfObjects=[]
listOfObjects.append(vessel)
listOfObjects.append(filling)
listOfObjects.append(wls)

for obj in listOfObjects:
    print(obj.name)
    obj.search_daughters(listOfObjects)
    print(environment.name)
    environment.search_daughters(listOfObjects)
    
# load settings into propagator
pro=Propagator(environment, listOfObjects)


##################
### ANGLE LOOP ###
##################

# use these zenith angles
zenith_range=np.linspace(0, np.pi/2, 90)

# save stuff from loop
Aeffs=[]
Aeff_errs=[]

# start loop
for zenith in zenith_range:
    print("zenith angle: {:.0f} [deg]".format(zenith*180/np.pi))
    
    # calc TIR analytically
    TIR, err=quad(TIR_Probability, 0,np.pi/2, args=(n_quartz, n_pfpe))

    # generate photons
    gen.rotatePlane = False
    gen.plane_zenith = zenith
    photons=gen.makeitso()
    
    # propagate photons
    count=0
    weights=0
    for ph in photons:
        ret=pro.loop(ph)
        weight=np.product(ph.weights)
        if ret==0:
            weights+=weight
            count+=1


    # uncertainties
    err_abs=np.sqrt(count)
    err_rel=err_abs/count

    Aplane=plane_x*plane_y
    Aeff=Aplane * weights / number
    print("Aeff temp: %.2f +/- %.2f cm^2"
        %(Aeff/100, Aeff/100*err_rel))

    Aeff*=all_factors*TIR
    print("Aeff: %.2f +/- %.2f cm^2"
        %(Aeff/100, Aeff/100*err_rel))
    Aeffs.append(Aeff/100)
    Aeff_errs.append(Aeff/100*err_rel)
    #break

# save data in csv file
angle = np.linspace(0, 90, 90)
dict = {"angle [deg]" : angle, "Aeff [cm2]" : Aeffs, "u(Aeff) [cm2]" : Aeff_errs}
df = pd.DataFrame(dict)
df.to_csv('efficiency_vs_zenith.csv')
