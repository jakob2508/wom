import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.integrate import quad
import pickle


def finvgaussian(x,params):
    lamb, mu = params
    y = np.heaviside(x-z/c*1E9,0) * np.sqrt(lamb/(2*np.pi*np.abs(x-z/c*1E9)**3))*np.exp(-lamb*((x-z/c*1E9)-mu)**2/(2*mu**2*(x-z/c*1E9)))
    y /= np.sum(y)
    return y

def loss(params):
    if params[0] < 0 or params[1] < 0:
        loss = 1E10
    else:
        loss = np.sqrt(np.sum((y-finvgaussian(x,params))**2))
    return loss

def integrant(xp, tau, x):
    return 1/xp**2 * np.exp(-(x-xp)/tau)

def fanalytic(x):
    tau = 1.6
    LoC= 0.38/299792458 * 1.47 * 1E9
    y = []
    for xx in x:
        y.append(np.heaviside(xx-z/c*1E9,0) * LoC/tau * quad(integrant, 1E-3, 20, args=(tau,xx))[0])
    y /= np.sum(y)
    return y

def func(x):
    y = np.heaviside(x-z/c*1E9,0)*z/c*1E9/x**2
    return y/np.sum(y)

# loading data
file = open("./pdf_wom_z=38cm.pkl","rb")
t, y = pickle.load(file)
file.close()
y2 = np.cumsum(y)

# inv gauss fit
x,z,c = t, 38/100, 299792458/1.47
res = minimize(loss, x0=(1,z/c*1E9), tol=1E-3, method="Nelder-Mead")
y_inv = finvgaussian(x,res.x)

# polinomial fit
fit = np.poly1d(np.polyfit(t, y, deg = 20))
y_poly = fit(x)

# analytical function
y_ana = fanalytic(x)

fig, ax = plt.subplots(1,1)
ax2 = plt.twinx(ax)
ax.step(t, y, where = "mid", color = "C0")
ax2.step(t,y2, where = "mid", color = "C1")
# ax.plot(x, y_inv, color = "C1", label = "inv. gaussian")
# ax.plot(x, y_poly, color = "red", label = "poly. 20th deg.")
# ax.plot(x, y_ana, color = "purple", label = "num. int.")
ax.set_xlabel("propagation time t [ns]")
ax.set_ylabel("pdf", color = "C0")
ax2.set_ylabel("cdf", color = "C1")
ax.set_xlim((0,20))
ax.set_ylim((0,0.035))
ax.tick_params(axis='y', colors='C0')
ax2.tick_params(axis='y', colors='C1')
ax2.set_ylim((0,1))
plt.tight_layout()
plt.show()
