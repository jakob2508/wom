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

def plot_3D_cone(ax, origin, phi, theta, length, opening_angle, color="blue"):

    phi = phi * np.pi/180
    theta = theta * np.pi/180
    # correct angles such that for (phi,theta) = (0,0) the cone is pointing along the positive x-axis direction
    phi = -phi
    # if theta < 0:
    #     theta -= 1*np.pi/4
    # elif theta > 0:
    theta -= np.pi/2

    # rotation matrix around axis y
    Ry = np.array([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta),0,np.cos(theta)]])
    # rotation matrix around axis z
    Rz = np.array([[np.cos(phi),-np.sin(phi),0],[np.sin(phi),np.cos(phi),0],[0,0,1]])

    # array of angles and heights
    t = np.linspace(0, 2 * np.pi, 100)
    u = np.linspace(0, length, 100)

    T, U = np.meshgrid(t,u)
    T, U = T, U

    X = U * np.tan(opening_angle) * np.cos(T)
    Y = U * np.tan(opening_angle) * np.sin(T) 
    Z = U

    # rotation and translation
    P = np.array([X,Y,Z])
    Py = P.T.dot(Ry).T 
    Pz = Py.T.dot(Rz).T
    Pt = Pz+origin[:, np.newaxis, np.newaxis]
    PP = Pt

    ax.plot_surface(PP[0], PP[1], PP[2],  rstride=4, cstride=4, color=color, linewidth=0, alpha=0.5)
    return ax