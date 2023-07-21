#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

#P = np.array([172, 160, 143, 124, 97, 86, 77, 75, 89, 125]) # en mW avant lentille
#U = np.array([2.27, 2.15, 1.93, 0.96, 0.74, 0.65, 0.56, 0.52, 0.547, 0.67]) # en volts

P = np.array([178,157,127,78,44,2, 43, 154, 140])
U = np.array([1.93,1.70,1.33,0.78,0.43,0.008, 0.40, 1.53, 1.39])

#plt.scatter(U,P)

def lin(x, a):
    return a*x

a,b = np.polyfit(U,P,1)
print(a,b)

a = opt.curve_fit(lin, U, P)[0][0]
print(a)

uu = np.linspace(0,3)

calib = {
        'xlabel': "Tension (mV)",
        'xdata': U,
        'xerr': 0.1,
        'xlims': (-0.5,2.5),
        'ylabel': "Puissance avant lentille (mW)",
        'ydata': P,
        'yerr': 5,
        'ylims': (-1,201),
        'func': lin,
        'fitlbl': "{:3.0f} mV/mw"
        #'params': a
        }


def plot_data(data):    
    fig, ax = plt.subplots()
    #fig = plt.figure()
    #ax = plt.axes((0.1,0.1,0.5,0.8))
    ax.errorbar(data['xdata'], data['ydata'], xerr=data['xerr'], yerr=data['yerr'], fmt='o')
    xx = np.linspace(data['xlims'][0], data['xlims'][1])
    a = opt.curve_fit(lin, U, P)[0]
    ax.plot(xx, data['func'](xx, *a), label=data['fitlbl'].format(*a))
    #ax.margins(0.3)
    #ax.use_sticky_edges = False
    ax.set_xlabel(data['xlabel'])
    ax.set_ylabel(data['ylabel'])
    ax.set_xlim(data['xlims'])
    ax.set_ylim(data['ylims'])

if False:
    plt.plot(uu, a*uu, label="{:3.0f} mV/mw".format(a), linestyle="--")
    plt.xlabel("Tension (mV)")
    plt.ylabel("Puissance avant lentille (mW)")

plot_data(calib)
plt.legend()
plt.show()

#fig, ax = plt.figure

#plt.savefig('calib photodiode.png', dpi=300)
