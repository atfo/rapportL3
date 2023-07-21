# -*- coding: utf-8 -*-
"""
Created on Wed May 31 15:43:29 2023

@author: alexandre.fouquet
"""

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

lmbd = 1.064 # um

z1 = np.array([2.5, 2.5+25.8-1.4, 105]) # distance au colimateur (cm)
wx1 = np.array([601.4, 646.9, 2239])/2  # waist horizontal (um)

z2 = np.array([9, 38, 95]) # distance au colimateur (cm)
wx2 = np.array([545.5, 684.2, 1978])/2  # waist horizontal (um)

z=np.hstack((z1,z2))
wx=np.hstack((wx1,wx2))

ordre=np.argsort(z)
z=z[ordre]; wx=wx[ordre]

w0 = 520 #um
zr = 28 #cm

def waist(z,w0,z0):
    zr=np.pi*w0**2/lmbd * 1e-4 # cm
    #print(zr)
    return w0*np.sqrt(1+(z-z0)**2 / zr**2)

popt, pcov = opt.curve_fit(waist, z, wx, bounds=([100,-500],[550,500]))

zz = np.linspace(-5,120)
plt.scatter(z,wx,label="waist mesuré")
plt.plot(zz,waist(zz,*popt), label="$w_0=$ {:.1f} $\\mathrm{{\\mu m}}$, $z_0=$ {:.1f} $\\mathrm{{cm}}$".format(*popt))
plt.legend()
plt.xlabel("z ($\mathrm{cm}$)")
plt.ylabel("waist ($\mathrm{\mu m}$)")
plt.title("Étude du mode en sortie de la fibre")
plt.figure()


plt.plot(z1,wx1,label="waist horiz")
plt.plot(z1,waist(z1,*popt))
plt.plot(z2,wx2,label="waist horiz")
plt.plot(z2,waist(z2,*popt))
plt.legend()
plt.ylabel("waist ($\mathrm{\mu m}$)")
plt.xlabel
plt.show()
