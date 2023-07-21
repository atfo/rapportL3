#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import scienceplots
import sys, argparse

import boyd

from boydf import alpha_np as balpha

import tikzplotlib
import matplotlib as mpl

plt.style.use('science')
mpl.rcParams['lines.linewidth'] = 3

def tpl_fix(obj):
    """
    workaround for matplotlib 3.6 renamed legend's _ncol to _ncols, which breaks tikzplotlib
    """
    if hasattr(obj, "_ncols"):
        obj._ncol = obj._ncols
    for child in obj.get_children():
        tpl_fix(child)

def savefig(name):
    tpl_fix(plt.gcf())
    tikzplotlib.save(name+".tex")

def idata(fname):
    return np.genfromtxt('./donnees/'+fname, delimiter=',', skip_header=1, usemask=True)

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
dcalib = idata('calib photodiode.csv')

#print(dcalib)

calib = {
        'name': 'calib photodiode',
        'xlabel': "Tension (V)",
        'xdata': dcalib[:,0], #[270, 271.8, 274.5, 277.4, 280, 283.3, 285.1, 286.6, 288.8, 290.6, 292.6, 297.3, 304.7, 307.8, 312.7, 324]
        'xerr': 0.01,
        'xlims': (0,1.6),
        'ylabel': "Puissance avant lentille (W)",
        'ydata': dcalib[:,1],
        'yerr': 1,
        'ylims': (0,5),
        'func': lin,
        'fitlbl': "{:3f} W/V"
        }

def quadra(x,a):
    return a* x**2

ddepp = idata('dep en basse puissance.csv')

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

def gprofile(z,w0,z0, lmbd=1.064):
    #lmbd = 1.064 # um
    zr=np.pi*w0**2/lmbd * 1e-4 # cm
    print('w0='+str(w0)+'zr='+str(zr))
    return w0*np.sqrt(1+(z-z0)**2 / zr**2)

# sortie de lentille
waist = {
        'name': 'waist faisceau incident',
        'xlabel': r"distance z à la lentille ($\mathrm{cm}$)",
        'xdata': 1.2+np.array([41, 26.5, 5, 3.5, 9.2, 3, 14, 27, 20]), # distance au bord de la lentille (cm), +1.2cm à cause de la profondeur caméra
        'xerr': 0.3,
        'xlims': (0,45),
        'ylabel': r"waist ($\mathrm{\mu m}$)",
        'ydata': np.array([3990, 2486, 200, 136, 602, 174, 1070, 2580, 1820])/2,  # waist horizontal (um)
        'yerr': np.array([200]+[80]*8)/2,
        'y2data': np.array([4100, 2740, 206, 133, 622, 184, 1230, 2780, 2040])/2,
        'ylims': (0,2000),
        'func': gprofile,
        'fitlbl': "$w_0=$ {:.1f} $\\mathrm{{\\mu m}}$, $z_0=$ {:.1f} $\\mathrm{{cm}}$",
        'datalbl': "waist mesuré",
        'nomask': True
        }

profil_vert = {
        'name': 'waist faisceau vert',
        'xlabel': r"distance z à la lentille ($\mathrm{cm}$)",
        'xdata': 1.2+np.array([39, 37.5, 34.1, 39, 34.5, 35, 35.5]), # distance au bord de la lentille (cm), +1.2cm à cause de la profondeur caméra
        'xerr': 0.1,
        'xlims': (0,45),
        'ylabel': r"waist ($\mathrm{\mu m}$)",
        'ydata': np.array([3543, 3505, 3000, 3638, 3262, 3211, 3300])/2,  # waist horizontal (um)
        'yerr': 50,
        'y2data': np.array([4100, 2740, 206, 133, 622, 184, 1230, 2780, 2040])/2,
        'ylims': (0,2000),
        'func': lambda z,w0,z0: gprofile(z,w0,z0,0.532),
        'fitlbl': "$w_0=$ {:.1f} $\\mathrm{{\\mu m}}$, $z_0=$ {:.1f} $\\mathrm{{cm}}$",
        'datalbl': "waist mesuré",
        'nomask': True
}

ddepTbP = idata('dep T basse P.csv')

pr = ddepTbP[:,0]
pv = ddepTbP[:,1]

pv = pv-10
neg = pv<0
pv[neg]=0

alpha = pv / (pr**2)

def boyd_fit(T,zr,lp,c):
    res = []
    for T in T:
        res.append(c*boyd.alpha(zr,T))
    return res



waist_shunte = {
        'name': 'waist faisceau laser',
        'xlabel': r"distance z au 2e miroir ($\mathrm{cm}$)",
        'xdata': 1.2+np.array([18.5, 7, -12, -24, 24.5]), # distance au bord de la lentille (cm), +1.2cm à cause de la profondeur caméra
        'xerr': 0.1,
        'xlims': (-30,80),
        'ylabel': r"waist ($\mathrm{\mu m}$)",
        'ydata': np.array([880, 845, 677, 785, 990])/2,  # waist horizontal (um)
        'yerr': 50,
        'ylims': (0,900),
        'func': gprofile,
        'fitlbl': "$w_0=$ {:.1f} $\\mathrm{{\\mu m}}$, $z_0=$ {:.1f} $\\mathrm{{cm}}$",
        'datalbl': "waist mesuré",
        'nomask': True
}

waist_incident_75 = {
        'name': 'waist faisceau incident (lentille 75 mm)',
        'xlabel': r"distance z à la lentille ($\mathrm{cm}$)",
        'xdata': 1.2+np.array([4.4,9.2,9.7,12.1,13.5, 12.8,11.8,7.9,8.5]), # distance au bord de la lentille (cm), +1.2cm à cause de la profondeur caméra
        'xerr': 0.1,
        'xlims': (0,16),
        'ylabel': r"waist ($\mathrm{\mu m}$)",
        'ydata': np.array([300,526,610,1120,1530, 1384, 1224, 293, 400])/2,  # waist horizontal (um)
        'yerr': 50,
        'ylims': (0,3500),
        'func': gprofile,
        'fitlbl': "$w_0=$ {:.1f} $\\mathrm{{\\mu m}}$, $z_0=$ {:.1f} $\\mathrm{{cm}}$",
        'datalbl': "waist mesuré",
        'nomask': True
}

waist_100 = {
        'name': 'waist faisceau incident (lentille 100 mm)',
        'xlabel': r"distance z à la lentille ($\mathrm{cm}$)",
        'xdata': 1.2+4+np.array([14,12,9,7.3,6,6.5,5.1,12.5,13.5,15.5,7]), # distance au bord de la lentille (cm), +1.2cm à cause de la profondeur caméra
        'xerr': 0.1,
        'xlims': (0,22),
        'ylabel': r"waist ($\mathrm{\mu m}$)",
        'ydata': np.array([1700,1225,930,505,234,276,90,1360,1566,1875,331])/2,  # waist horizontal (um)
        'yerr': 50,
        'ylims': (0,2000),
        'func': gprofile,
        'fitlbl': "$w_0=$ {:.1f} $\\mathrm{{\\mu m}}$, $z_0=$ {:.1f} $\\mathrm{{cm}}$",
        'datalbl': "waist mesuré",
        'nomask': True
}

waist_150 = {
        'name': 'waist faisceau incident (lentille 150 mm)',
        'xlabel': r"distance z à la lentille ($\mathrm{cm}$)",
        'xdata': 1.2+np.array([6.5,10.8,14.5,17]), # distance au bord de la lentille (cm), +1.2cm à cause de la profondeur caméra
        'xerr': 0.1,
        'xlims': (0,22),
        'ylabel': r"waist ($\mathrm{\mu m}$)",
        'ydata': np.array([286,312,644,875])/2,  # waist horizontal (um)
        'yerr': 50,
        'ylims': (0,2000),
        'func': gprofile,
        'fitlbl': "$w_0=$ {:.1f} $\\mathrm{{\\mu m}}$, $z_0=$ {:.1f} $\\mathrm{{cm}}$",
        'datalbl': "waist mesuré",
        'nomask': True
}

def plot_data(data, p=None):    
    fig, ax = plt.subplots()#(2,1)
    #ax = axs[1]
    #fig = plt.figure()
    #ax = plt.axes((0.1,0.1,0.5,0.8))
    x = data['xdata']
    y = data['ydata']
    ax.errorbar(data['xdata'], data['ydata'], xerr=data['xerr'], yerr=data['yerr'], fmt='.', markersize=0.2, label=data['datalbl'])
    if 'func' in data and p==None:    
        xx = np.linspace(data['xlims'][0], data['xlims'][1])
        if not ('nomask' in data):
            x,y = x[~x.mask], y[~x.mask]
        if 'p0' in data:
            a = opt.curve_fit(data['func'], x, y, p0=data['p0'])[0]
        else:
            a = opt.curve_fit(data['func'], x, y)[0]
        ax.plot(xx, data['func'](xx, *a), label=data['fitlbl'].format(*a))
        print(*a)
    if p:
        xx = np.linspace(data['xlims'][0], data['xlims'][1], 1000)
        ax.plot(xx, data['func'](xx, *p), label=data['fitlbl'].format(*p))
        #axs[0].plot(xx, data['func'](xx,*p))
    #ax.margins(0.3)
    #ax.use_sticky_edges = False
    ax.set_xlabel(data['xlabel'])
    ax.set_ylabel(data['ylabel'])
    ax.set_xlim(data['xlims'])
    ax.set_ylim(data['ylims'])
    ax.legend()
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', action='store_true')    
    if parser.parse_args().s:    
        #plt.savefig(data['name']+".pdf", dpi=300)
        savefig(data['name'])

if False:
    plt.plot(uu, a*uu, label="{:3.0f} mV/mw".format(a), linestyle="--")
    plt.xlabel("Tension (mV)")
    plt.ylabel("Puissance avant lentille (mW)")

#plot_data(calib)
#plot_data(depp)
#plot_data(waist)
#plot_data(depTbP)
#plot_data(depTbP, p=[1.2,6.8741,0.143])
#plot_data(waist_100)


waist_test = {
        'name': 'waist faisceau incident (lentille 100mm)',
        'xlabel': r"distance z à la lentille ($\mathrm{cm}$)",
        'xdata': 1.2+np.array([7.3,6,4.5,7.5,11.2,12.7]), # distance au bord de la lentille (cm), +1.2cm à cause de la profondeur caméra
        'xerr': 0.1,
        'xlims': (0,22),
        'ylabel': r"waist ($\mathrm{\mu m}$)",
        'ydata': np.array( [681,400,128,671,1360,1665])/2,  # waist horizontal (um)
        'yerr': 50,
        'ylims': (0,2000),
        'func': gprofile,
        'fitlbl': "$w_0=$ {:.1f} $\\mathrm{{\\mu m}}$, $z_0=$ {:.1f} $\\mathrm{{cm}}$",
        'datalbl': "waist mesuré",
        'nomask': True
}

#plot_data(waist_150)

#plot_data(profil_vert)


depTbP = {
    'name': 'Efficacité de conversion à basse puissance en fonction de la température',
    'xlabel': "Température (°C)",
    'xdata': ddepTbP[:,2],
    'xlims': (77.8,90.2),
    'ylabel': r'$\alpha$ $(W/W^2)$',
    'ydata': alpha,
    'ylims': (-0.001,0.04),
    'xerr': 0.001,
    'yerr': 0.0001,
    'func': lambda T,zr,lp,c: balpha(zr,T,lp,c),
    'p0': [1.2,6.8741,0.143],
    'bounds': ([0.5,]),
    'fitlbl': 'Ajustement'#' à la théorie de Boyd-Kleinman (zr = {}, lp = {}, c = {})'
}

d75mf = idata('f75malfoc.csv')

print(d75mf[:,0]+8)

depT75mf = {
    'name': 'Efficacité de conversion à basse puissance en fonction de la température (f = 75 mm)',
    'xlabel': "Température (°C)",
    'xdata': d75mf[:,0]+8,
    'xlims': (77.8,90.2),
    'ylabel': r'$\alpha$ $(W/W^2)$',
    'ydata': d75mf[:,1]*1e-2,
    'ylims': (-0.001,0.04),
    'xerr': 0.001,
    'yerr': 0.0001,
    'func': lambda T,zr,lp,c: balpha(zr,T,lp,c),
    'p0': [1.2,6.8741,0.143],
    'bounds': ([0.5,]),
    'fitlbl': 'Ajustement'#' à la théorie de Boyd-Kleinman (zr = {}, lp = {}, c = {})'
}

plot_data(waist)
plot_data(profil_vert)
plt.figure()
plt.scatter(waist['xdata'],waist['ydata'], label='fondamental')
plt.scatter(profil_vert['xdata'],profil_vert['ydata'], label='seconde harmonique')
plt.legend()
plt.xlabel(waist['xlabel'])
plt.ylabel(waist['ylabel'])

#plot_data(depTbP)
#plot_data(depTbP, p=[2.48,6.862,0.097*(17/14)**2])
#plot_data(depT75mf, p=[1.7,6.862,0.097*(17/14)**2])
savefig("comparaison-profils")
plt.show()


