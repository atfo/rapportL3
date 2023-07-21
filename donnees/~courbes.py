#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import scienceplots

from .. import boyd

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

#dcalib = np.genfromtxt('calib photodiode 9 juin.csv', delimiter=',', skip_header=1, usemask=True)
#dcalib = np.genfromtxt('calib glan.csv', delimiter=',', skip_header=1, usemask=True)
dcalib = np.genfromtxt('calib photodiode 15 juin.csv', delimiter=',', skip_header=1, usemask=True)

#print(dcalib)

calib = {
        'name': 'calib photodiode',
        'xlabel': "Tension (V)",
        'xdata': dcalib[:,0], #[270, 271.8, 274.5, 277.4, 280, 283.3, 285.1, 286.6, 288.8, 290.6, 292.6, 297.3, 304.7, 307.8, 312.7, 324]
        'xerr': 0.01,
        'xlims': (0,1.5),
        'ylabel': "Puissance avant lentille (mW)",
        'ydata': dcalib[:,1],
        'yerr': 1,
        'ylims': (0,150),
        'func': lin,
        'fitlbl': "{:3.0f} mW/V"
        }

def quadra(x,a):
    return a* x**2

ddepp = np.genfromtxt('dep en basse puissance.csv', delimiter=',', skip_header=1, usemask=True)

conv = 233
depp = {
        'name': 'conversion basse puissance 82.5 C',
        'xlabel': "Puissance en entrée (mW)",
        'xdata': ddepp[:,0]*conv, # conversion des V à la diode en mW en entrée
        'xerr': 0.01*conv,
        'xlims': (0,250),
        'ylabel': "Puissance de vert ($\\mathrm{{\\mu W}}$)",
        'ydata': ddepp[:,1],
        'yerr': 0.5,
        'ylims': (0, 2000),
        'func': quadra,
        'fitlbl': "$\\alpha = \\frac{{P_2}}{{P_1^2}} = {:.3} \mathrm{{W/W^2}}$"
}

def gprofile(z,w0,z0):
    lmbd = 1.064 # um
    zr=np.pi*w0**2/lmbd * 1e-4 # cm
    #print(zr)
    return w0*np.sqrt(1+(z-z0)**2 / zr**2)

# sortie de lentille
waist = {
        'name': 'waist faisceau incident',
        'xlabel': r"distance z à la lentille ($\mathrm{cm}$)",
        'xdata': 1.2+np.array([41, 26.5, 5, 3.5, 9.2, 3, 14, 27, 20]), # distance au bord de la lentille (cm), +1.2cm à cause de la profondeur caméra
        'xerr': 0.1,
        'xlims': (0,45),
        'ylabel': r"waist ($\mathrm{\mu m}$)",
        'ydata': np.array([3990, 2486, 200, 136, 602, 174, 1070, 2580, 1820])/2,  # waist horizontal (um)
        'yerr': np.array([200]+[20]*8)/2,
        'y2data': np.array([4100, 2740, 206, 133, 622, 184, 1230, 2780, 2040])/2,
        'ylims': (0,4000),
        'func': gprofile,
        'fitlbl': "$w_0=$ {:.1f} $\\mathrm{{\\mu m}}$, $z_0=$ {:.1f} $\\mathrm{{cm}}$",
        'datalbl': "waist mesuré",
        'nomask': True
        }
ddepTbP = np.genfromtxt('dep T basse P.csv', delimiter=',', skip_header=1, usemask=True)

pr = ddepTbP[:,0]
pv = ddepTbP[:,1]

pv = pv-13
neg = pv<0
pv[neg]=0

alpha = pv / (pr**2)

depTbP = {
    'name': 'Efficacité de conversion à basse puissance en fonction de la température',
    'xlabel': "Température (°C)",
    'xdata': ddepTbP[:,2],
    'xlims': (77.8,90.2),
    'ylabel': r'$\alpha$ $(W/W^2)$',
    'ydata': alpha,
    'ylims': (0,0.04),
    'xerr': 0.001,
    'yerr': 0.0001,
    'func': lambda T,zr,c: c*boyd.alpha(zr,T),
    'fitlbl': 'Ajustement à la théorie de Boyd-Kleinman'
}

def plot_data(data):    
    fig, ax = plt.subplots()
    #fig = plt.figure()
    #ax = plt.axes((0.1,0.1,0.5,0.8))
    x = data['xdata']
    y = data['ydata']
    ax.errorbar(data['xdata'], data['ydata'], xerr=data['xerr'], yerr=data['yerr'], fmt='.', markersize=0.2)
    if 'func' in data:    
        xx = np.linspace(data['xlims'][0], data['xlims'][1])
        if 'nomask' in data:
            a = opt.curve_fit(data['func'], x, y)[0]
        else:
            a = opt.curve_fit(data['func'], x[~x.mask], y[~x.mask])[0]
        ax.plot(xx, data['func'](xx, *a), label=data['fitlbl'].format(*a))
    #ax.margins(0.3)
    #ax.use_sticky_edges = False
    ax.set_xlabel(data['xlabel'])
    ax.set_ylabel(data['ylabel'])
    ax.set_xlim(data['xlims'])
    ax.set_ylim(data['ylims'])
    ax.legend()
    #plt.savefig(data['name']+".pdf", dpi=600)

if False:
    plt.plot(uu, a*uu, label="{:3.0f} mV/mw".format(a), linestyle="--")
    plt.xlabel("Tension (mV)")
    plt.ylabel("Puissance avant lentille (mW)")

plot_data(calib)
plot_data(depp)
plot_data(waist)
plot_data(depTbP)
plt.show()


