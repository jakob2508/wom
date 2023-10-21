import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.optimize import minimize

# We fit the efficiency data with a linear function f(z) = a*z + b.
# However, there are constraints to the parameters a and b because the sum of both PMTs is always 1 = 100%.
# We know that the sum f_PMT1(z) + f_PMT2(z) = 1 and that the function is symmetric f_PMT(z) = f_PMT(height-z)
# That we can write f_PMT1(0) = b
# And a*height + b = f_PMT1(height) = 1 - f_PMT2(height) = 1 - f_PMT1(height-heigth) = 1 - b
# Thus we can replace b = (1 - a*height)/2 and effectively fit one parameter!
def linear_func(x,args):
    a = args
    return a*x+((1-a*height)/2)

def fit_func(arg):
    chi2 = np.sqrt(np.sum((y-linear_func(x,arg))**2/y_err**2))
    return chi2

# load the data
file = open("./files/sim_eff.pkl","rb")
z_range, eff1, u_eff1, eff2, u_eff2 = pickle.load(file)

# efficiency in %
eff1, u_eff1 = eff1*100, u_eff1*100
eff2, u_eff2 = eff2*100, u_eff2*100
efftot, u_efftot = eff1+eff2, np.sqrt((eff2*u_eff1)**2+(eff1*u_eff2))

# relative efficiency = PMT selection probability
# losses due to abosorption already included in the eff. area
reff1 = eff1/efftot
reff2 = eff2/efftot
refftot = reff1+reff2

u_reff1 = np.sqrt((u_eff1/efftot)**2+(u_efftot*eff1/efftot**2)**2)
u_reff2 = np.sqrt((u_eff2/efftot)**2+(u_efftot*eff2/efftot**2)**2)

height = 76 # inner cylinder height

# fit the renormalized efficiency data with a linear function
# for the upper PMT
x, y, y_err = z_range, reff1, u_reff1
res1 = minimize(fit_func, x0 = (1))
# and the lower PMT
x, y, y_err = z_range, reff2, u_reff2
res2 = minimize(fit_func, x0 = (1))

# plot of the effciency as a function of the position z
fig, ax = plt.subplots(1,2, figsize = (10,4))
ax[0].plot(z_range, eff1, color = "grey", linestyle = ":", label = r"$\epsilon_{PMT,up}$")
ax[0].plot(z_range, eff2, color = "black", linestyle = "--", label = r"$\epsilon_{PMT,down}$")
ax[0].plot(z_range, efftot, color = "red", label = r"$\epsilon_{tot}$")
ax[0].set_xlabel("z position [cm]")
ax[0].set_ylabel(r"efficiency $\epsilon$ [%]")
ax[0].set_ylim((0,30))
ax[0].legend()

# plot of the PMT selection probability as a function of the position z
ax[1].errorbar(z_range, reff1*100, yerr = u_reff1*100, color = "grey", marker = "o", markersize = 3, linestyle = "none", capsize = 3, label = "PMT,up")
ax[1].plot(z_range, linear_func(z_range, res1.x[0])*100, color = "grey", label = "fit: a={:.3f} 1/m".format(res1.x[0]*100))
ax[1].errorbar(z_range, reff2*100, yerr = u_reff2*100, color = "black", marker = "s", markersize = 3, linestyle = "none", capsize = 3, label = "PMT,down")
ax[1].plot(z_range, linear_func(z_range, res2.x[0])*100, color = "black", label = "fit: a={:.3f} 1/m".format(res2.x[0]*100))
ax[1].plot(z_range, refftot*100, color = "red", label = "100%")
ax[1].set_xlabel("z position [cm]")
ax[1].set_ylabel(r"PMT hit probability [%]")
ax[1].set_ylim((0,105))
ax[1].legend()

plt.tight_layout()
plt.savefig("./plots/loss_vs_distance.png")
plt.show()

# We take the average slope of both efficiency fits for the upper and lower PMT.
av_a = np.mean([res1.x[0], -res2.x[0]])

print("Fit parameter PMT_up/PMT1, y=a*x+(1-ah)/2: a={:.6f}".format(res1.x[0]*100))
print("Fit parameter PMT_down/PMT0, y=a*x+(1-ah)/2: a={:.6f}".format(res2.x[0]*100))
print("Fit parameter combined result, y=a*x+(1-ah)/2: a={:.6f}".format(av_a*100))
