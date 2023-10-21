import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

# load data
file = pd.read_csv("efficiency_vs_zenith.csv", sep = ',' , names=["x","y","z"], skiprows=1)
angle = file["x"].values # angle
Aeffs = file["y"].values # effective area
Aeff_errs = file["z"].values # error on effective area

Aeff_errs = np.where(Aeffs == 0, 1E-3, Aeff_errs) # removing nans from Aeff_errs

# zenith angle in radians of plane's normal vector relative to the cylinder's symmetry axis
zenith=np.linspace(np.pi/2, 0, 90)

# Mirrorign negative angles: f(-zenith) = f(zenith)
# cos(zenith), x value
coss = np.stack((np.cos(zenith), -np.cos(zenith))).flatten()
# effective area, integral normalized to 1
pdfs = np.stack((Aeffs/np.sum(Aeffs), Aeffs/np.sum(Aeffs))).flatten()
pdfs = pdfs * 1/np.max(pdfs) * 2/np.pi # 2/pi normalization
# uncertainty on effective area
pdfs_err = np.stack((Aeff_errs/np.sum(Aeffs), Aeff_errs/np.sum(Aeffs))).flatten()

# sort data from -pi/2 to pi/2
sort = np.argsort(coss)
coss = coss[sort]
pdfs = pdfs[sort]
pdfs_err = pdfs_err[sort]

x = np.linspace(-1,1,len(coss))
# sin(theta) as reference
sins = np.sin(np.arccos(x))*2/np.pi

# fit polynomial to effective area as a function of cos(theta), returns coefficients of a 10th degree polynimial
coef = np.polyfit(coss, pdfs, w = pdfs_err, deg = 10)
fit = np.poly1d(np.polyfit(coss, pdfs, w = pdfs_err, deg = 10))

# plotting
fig, ax = plt.subplots(1,2, figsize = (10,4))
fig.subplots_adjust(wspace = 0.3)
ax[0].errorbar(zenith*180/np.pi, Aeffs, yerr = Aeff_errs, fmt = "None", capsize = 2, label = "MC")
ax[0].plot(zenith*180/np.pi, np.max(Aeffs)*np.sin(zenith), label = r"$sin(\theta)$", color = "C1")
ax[0].set_xlabel(r"$\theta$")
ax[0].set_ylabel(r"effective area $A_{eff} [cm^2]$")
ax[0].legend()
ax[1].errorbar(coss, pdfs, yerr = pdfs_err, fmt = "None", capsize = 2, label = "MC")
ax[1].plot(x,sins, label = r"$sin(\theta)$", color = "C1")
ax[1].plot(x, fit(x), label = "fit", color = "C2")
ax[1].set_xlabel(r"$cos(\theta)$")
ax[1].set_ylabel(r"angular sensitivity PDF")
ax[1].legend()

plt.tight_layout()
plt.savefig("angular_sensitivity_wom.png")
plt.show()


# Angular sensitivity file contains the maximum value of the angular sensitivity curve, as well as the polinomial coefficients of 10th order to form p0, p1, p2, ... where p0*x**0 + p1*x**1 + p2*x**2 ...
if 1:
    filename = "as_140.0"
    file = open(filename, 'w')
    sys.stdout = file
    print("{:f}".format(pdfs.max()))
    for i in range(len(coef)):
        if np.abs(coef[i]) < 1E-5:
            print("0")
        else:
            print("{:f}".format(coef[i]))
    file.close()
sys.stdout = sys.__stdout__
print('Finished!')