import numpy as np
import matplotlib.pyplot as plt

def plot_3D_cylinder(ax, rho, height, x0 = 0, y0 = 0, z0 = 0, color='k'):
    rho *= 0.99
    x = np.linspace(x0-rho, x0+rho, 100)
    z = np.linspace(z0-height/2, z0+height/2, 100)
    X, Z = np.meshgrid(x, z)

    Y = np.sqrt(rho**2 - (X-x0)**2) + y0 # Pythagorean theorem

    ax.plot_surface(X, Y, Z, linewidth=0, color=color, alpha = 0.5)
    ax.plot_surface(X, (2*y0-Y), Z, linewidth=0, color=color, alpha = 0.5)
    return ax

def plot_3D_sphere(ax, radius, x0, y0, z0, color='red'):
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = x0 - radius * np.outer(np.cos(u), np.sin(v))
    y = y0 - radius * np.outer(np.sin(u), np.sin(v))
    z = z0 - radius * np.outer(np.ones(np.size(u)), np.cos(v))

    ax.plot_surface(x, y, z,  rstride=4, cstride=4, color=color, linewidth=0, alpha=0.5)
    return ax