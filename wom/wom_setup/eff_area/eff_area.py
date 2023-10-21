import sys
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def folding_function(x):
    # 1/lambda**2 folding function
    f = 1/x**2
    norm = 1/245-1/425 # normalization = 1 in wavelengths interval [245, 425] nm
    return f/norm

# Convolution of two functions
def convolution(x, func1, func2):
    conv = sum([func2[(x-i)]*func1[i] for i in range(len(func1))])
    return conv


debug1 = 0 # eff area vs wavelength for wom in air
debug2 = 0 # eff area vs wavelength for wom in ice, folded with Cherenkov spectrum

file = pd.read_csv("eff_area_wom_air.csv", sep = ',' , names=["wvl","ea"])
wavelength = file["wvl"]
eff_area_air = file["ea"]

# The data in the csv file is for the WOM embedded in air. In ice however a maximum effective area of 18 cm^2 is reached.
max_eff_area_ice = 18 # 18 cm2 maximum eff area in ice
medium_scaling = max_eff_area_ice/np.max(eff_area_air)

# Only the inner tube is sensitive, therefore we adapt the height of the tube in PPC but rescale the eff. area by the ratio of the two heights. If the tube is shortened but measures as many photons we have to increase the efficiency/ effective area.
h_vessel = 130 # 130 cm vessel height
h_inner = 76 # 76 cm inner tube height
geometry_scaling = h_vessel/h_inner

wavelength_range = np.arange(245, 425, 5)

if debug1:
    fig,ax = plt.subplots(1,1)
    ax.plot(wavelength, eff_area_air, '.')
    ax.set_xlabel("wavelength [nm]")
    ax.set_ylabel(r"WOM effective area in air [cm$^2$]")
    plt.show()

# interpolation of effective area
f_eff_area = interp1d(wavelength, eff_area_air, kind = "quadratic")

# rescale effective area
eff_area_ice = f_eff_area(wavelength_range) * medium_scaling * geometry_scaling

# just for plotting, convolute rescaled eff area with Cherenkov spectrum
conv = convolution(np.arange(len(wavelength_range)), eff_area_ice, folding_function(wavelength_range))

if debug2:
    fig,ax = plt.subplots(1,1)
    ax.plot(wavelength_range, eff_area_ice, label = "rescaled eff. area")
    # scale Cherenkov folded eff area to arbitrary value to match scale
    ax.plot(wavelength_range, eff_area_ice*folding_function(wavelength_range)/folding_function(wavelength_range[0]), label = "Cherenkov weighted eff. area [AU]")
    # scale Cherenkov spectrum to arbitrary value to match scale
    ax.plot(wavelength_range, folding_function(wavelength_range)*30/np.max(folding_function(wavelength_range)), color = "grey", linestyle="--", label = "Cherenkov spectrum [AU]")
    ax.set_xlabel("wavelength [nm]", fontsize = 14)
    ax.set_ylabel(r"WOM effective area in ice [cm$^2$]", fontsize = 14)
    ax.tick_params(labelsize=14)
    plt.legend()
    plt.tight_layout()
    plt.savefig("eff_area_wom.png")
    plt.show()


if 1:
    filename = "om.wv_140.0"
    file = open(filename, 'w')
    sys.stdout = file
    for i in range(len(wavelength_range)):
        print("{:.0f} {:f}".format(wavelength_range[i],eff_area_ice[i]))

    file.close()
sys.stdout = sys.__stdout__
print('Finished!')
