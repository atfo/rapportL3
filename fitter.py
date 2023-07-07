#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import scienceplots
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from tqdm import tqdm
import sys, argparse

import boyd

def idata(fname):
    return np.genfromtxt('./donnees/'+fname, delimiter=',', skip_header=1, usemask=True)

ddepTbP = idata('dep T basse P.csv')

pr = ddepTbP[:,0]
pv = ddepTbP[:,1]

pv = pv-13
neg = pv<0
pv[neg]=0
T = ddepTbP[:,2]
alpha = pv / (pr**2)

def boyd_fit(T,zr,lp,c):
    return c*boyd.alpha(zr,T,lp)

depTbP = {
    'name': 'Efficacité de conversion à basse puissance en fonction de la température',
    'xlabel': "Température (°C)",
    'xdata': T,
    'xlims': (77.8,90.2),
    'ylabel': r'$\alpha$ $(W/W^2)$',
    'ydata': alpha,
    'ylims': (-0.001,0.04),
    'xerr': 0.001,
    'yerr': 0.0001,
    'func': lambda T,zr,lp,c: c*boyd.alpha(zr,T,lp),
    'p0': [1.2,6.9,0.07],
    'fitlbl': 'Ajustement à la théorie de Boyd-Kleinman'
}

def sqe(y,func,x,*p):
    #rs = (y-func(x,*p))**2
    return np.mean(np.square(y-func(x,*p)))

zR_arr = np.linspace(0.01, 1.5, 201)
T_arr = np.linspace(50, 90, 201)
lp_arr = np.linspace(6.8735,6.8745,31)
c_arr = np.linspace(0.122, 0.124, 5)

def fit(y,func,x,a,b):
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', action='store_true')    
    if parser.parse_args().u:
        na, nb = a.size, b.size
        im = np.empty((na, nb))
        for i in tqdm(range(na)):
            for j in range(nb):
                im[i,j] = sqe(y, func, x, 1.2, a[i], b[j])
        np.savez("fitter", im=im)
    else:
        im = np.load("fitter.npz")['im']
    plt.imshow(im, extent=(np.amin(b), np.amax(b), np.amin(a), np.amax(a)), aspect="auto", cmap=cm.inferno) #, norm=LogNorm())


#parser = argparse.ArgumentParser()
#parser.add_argument('-u', action='store_true')    
#if parser.parse_args().u:
#    im = np.empty((41, 21))
#    for i in tqdm(range(41)):
#        for j in range(21):
#            im[i,j] = sqe(alpha, boyd_fit, T, 1.2, lp_arr[i], c_arr[j])
#    np.savez("fitter", im=im)
#else:
#    im = np.load("fitter.npz")['im']
#
#plt.imshow(im, extent=(np.amin(c_arr), np.amax(c_arr), np.amin(lp_arr), np.amax(lp_arr)), aspect="auto", cmap=cm.inferno) #, norm=LogNorm())



fit(alpha, boyd_fit, T, lp_arr, c_arr)
print(sqe(alpha, boyd_fit, T, 1.2, 6.874, 0.123))


plt.xlabel('Coefficient')
plt.ylabel(r'$\Lambda$ ($\mu m$)')
plt.colorbar()
plt.show()

