# -*- coding: utf-8 -*-
"""
Created on Wed May 31 15:43:29 2023

@author: alexandre.fouquet
"""

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

lmbd = 1.064 # um

# sortie laser
#z1 = np.array([30, 13, 60, 92]) # distance au colimateur (cm)
#wx1 = np.array([905, 840, 990, 1260])/2  # waist horizontal (um)

# sortie telescope avant fibre
#z1 = np.array([92, 43.5, 56.5]) # distance au colimateur (cm)
#wx1 = np.array([1030, 650, 645])/2  # waist horizontal (um)

#sortie fibre
#z1 = np.array([7, 15, 88, 56]) # distance au colimateur (cm)
#wx1 = np.array([500, 536, 2400, 1730])/2  # waist horizontal (um)

#z1 = np.array([62,17,40]) # distance au colimateur (cm)
#wx1 = np.array([3600,360,2000])/2  # waist horizontal (um)

# sortie de lentille
z1 = np.array([41, 26.5, 5, 3.5, 9.2, 3, 14, 27, 20]) # distance au bord de la lentille (cm)
wx1 = np.array([3990, 2486, 200, 136, 602, 174, 1070, 2580, 1820])/2  # waist horizontal (um)
dwx1 = np.array([200, 5, 5, 5,5,5])/2
wy1 = np.array([4100, 2740, 206, 133, 622, 184, 1230, 2780, 2040])/2

z2 = np.array([]) # distance au colimateur (cm)
wx2 = np.array([])/2  # waist horizontal (um)

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

popt, pcov = opt.curve_fit(waist, z, wx)

zz=np.linspace(-5,50)

plt.scatter(z,wx,label="waist horiz")
plt.plot(zz,waist(zz,*popt), label="$w_0=$ {}, $z_0=$ {}".format(*popt))
plt.legend()
plt.figure()

#
#plt.plot(z1,wx1,label="waist horiz")
#plt.plot(z1,waist(z1,*popt))
#plt.plot(z2,wx2,label="waist horiz")
#plt.plot(z2,waist(z2,*popt))
#plt.legend()
#plt.show()