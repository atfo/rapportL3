# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 13:55:05 2023

@author: Alexandre
"""

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

#P = np.array([172, 160, 143, 124, 97, 86, 77, 75, 89, 125]) # en mW avant lentille
#U = np.array([2.27, 2.15, 1.93, 0.96, 0.74, 0.65, 0.56, 0.52, 0.547, 0.67]) # en volts

P = np.array([178,157,127,78,44,2, 43, 154, 140])
U = np.array([1.93,1.70,1.33,0.78,0.43,0.008, 0.40, 1.53, 1.39])

plt.scatter(U,P)

def lin(x, a):
    return a*x

a,b = np.polyfit(U,P,1)
print(a,b)

a = opt.curve_fit(lin, U, P)[0]
print(a)

uu = np.linspace(0,3)

plt.plot(uu, a*uu)
plt.xlabel("Tension (mV)")
plt.ylabel("Puissance avant lentille (mW)")