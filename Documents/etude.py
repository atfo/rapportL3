#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy.integrate as integr
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from tqdm import tqdm
import sys, argparse
import matplotlib as mpl
import tikzplotlib
import scienceplots
import scipy.stats as st

from scipy.optimize import curve_fit

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

L=2 #cm
lmbd1 = 1.064 #um
lmbd2 = 0.532

cexp, T0 = 1.5e-5, 20 # coeff dilat thermique


bb = np.linspace(-20,20,1000)
TT = np.linspace(50,150,500)

def idata(fname, mask=True):
    return np.genfromtxt('./'+fname, delimiter=',', skip_header=1, usemask=mask)

def n(T,lmbd):
    a = np.array([5.756, 0.0983, 0.2020, 189.32, 12.52, 1.32e-2])#.reshape(6,1)
    b = np.array([2.860e-6, 4.7e-8, 6.113e-8, 1.516e-4])#.reshape(4,1)
    f = (T-24.5)*(T+570.82)
    p = a[:4] + f*b
    inds = np.sqrt( p[0] + p[1]/(lmbd**2-p[2]**2) + p[3]/(lmbd**2-a[4]**2) - a[5]*lmbd**2)
    return inds

def n1(T):
    return n(T, lmbd1)

def n2(T):
    return n(T, lmbd2)

def bT(T,lp,zr):
    lpT = lp*(1+cexp*(T-20))
    return zr*2*np.pi* (1/lpT - (n2(T)-n1(T))/lmbd2 ) * 1e4

def Tdeb(b,lp,zr):
    TT = np.linspace(50,150,500)
    arr = np.array([bT(T,lp,zr) for T in TT])
    return TT[arr<b][0] # T(b) décr

def dkeff(T,lp):
    ''' lp en um, dkeff en cm-1 '''
    lpT = lp*(1+cexp*(T-20))
    return 2*np.pi* (1/lpT - (n2(T)-n1(T))/lmbd2 ) * 1e4

def lpopt(T):
    return 2*np.pi / (1.6e-4 + 2*np.pi*(n2(T)-n1(T))/lmbd2) / (1+cexp*(T-T0))

def h(a,b):
    def f(t):
        return np.exp(1j*b*t)/(1+1j*t)
    N = 1000
    x_array = np.linspace(-a, a, N).reshape(N,1)
    f_array = np.array([f(x) for x in x_array])
    f_array = f_array.reshape(N, f_array.size // N)
    #print(f_array.shape)
    integral = np.trapz(f_array, x_array, axis=0)
    return 1 / (4*a) * np.square(np.abs(integral))

#print(h(2.8,0.56))

def h_shift(a,b):
    ''' Si focalisé au bord du cristal comme on le craint '''
    def f(t):
        return np.exp(1j*b*t)/(1+1j*t)
    N = 1000
    x_array = np.linspace(-2*a, 0, N).reshape(N,1)
    f_array = np.array([f(x) for x in x_array])
    f_array = f_array.reshape(N, f_array.size // N)
    #print(f_array.shape)
    integral = np.trapz(f_array, x_array, axis=0)
    return 1 / (4*a) * np.square(np.abs(integral))

def hmax(a):
    return np.max(h(a,bb))

def bopt(a):
    return bb[np.argmax(h(a,bb))]

def fwhmf(x,y):
    #y = np.array(y)
    inds = y>np.max(y)/2
    return (x[inds][0],x[inds][-1])

def fwhmb(a):
    return fwhmf(bb, h(a,bb))

def caract(zr,lp):
    a = L/2/zr
    bo = bopt(a)
    hm = hmax(a)
    bm,bp = fwhmb(a)
    db = bp-bm
    dkeff_opt = - bo/zr # cm-1
    Topt = Tdeb(bo,lp,zr)
    dT = Tdeb(bp,lp,zr)-Tdeb(bm,lp,zr)
    lp_opt = 2*np.pi / (2*np.pi*(n2(Topt)-n1(Topt))/lmbd2 - dkeff_opt) / (1+cexp*(Topt-T0)) # à T0
    c= {
        'a' : a,
        'bo' : bo,
        'hm' : hm,
        'bm,bp' : (bm,bp),
        'db' : bp-bm,
        'dkeff_opt' : - bo/zr, # cm-1
        'Topt' : Topt,
        'dT' : dT,
        'lp' : lp_opt
        #dT2 = - 2*np.pi * zr /lmbd2 * 1e4 * (n2(100)-n1(100)-n2(80)+n1(80))/20
    }
    return c

#caract(2.4,6.9)


if __name__ == "__main__":
    if False:
        ddepp = idata('dep en basse puissance.csv',mask=False)
        conv = 233
        puissances = ddepp[:,0]*conv
        vert = ddepp[:,1] * 1e-3
        pe = 0.01 * conv
        ve = 2e-2
        print(vert)
        def quadra(x, a):
            return a * x ** 2
        popt,pcov = curve_fit(quadra,puissances,vert,p0=[0.02])
        with plt.style.context(['science','scatter']):
            plt.errorbar(puissances,vert,xerr=pe,yerr=ve,marker=',')
            p_arr = np.linspace(0,350,500)
            plt.plot(p_arr, quadra(p_arr,*popt), 'g--', label=r'$\alpha={}$'.format(popt[0]))
        plt.legend()
        savefig("quadra")
        plt.show()
    #exit(0)

    if False:
        for zr in [0.5,1,1.5,2,2.5]:
            T_arr = np.linspace(50,90,500)
            #plt.plot(T_arr, [h(1/zr, bT(T,6.9,zr)) for T in T_arr])
            #plt.show()
            dT = (fwhmf(T_arr, [h(1/zr, bT(T,6.9,zr))[0] for T in T_arr]))
            print(zr,dT, dT[1]-dT[0])
            #plt.show()

    def interval(data):
        return st.t.interval(alpha=0.2, df=len(data) - 1, loc=np.mean(data), scale=st.sem(data))

    if True:
        depT = idata('dep T basse P.csv', mask=False)
        dataT = depT[:,3]
        dataa = depT[:,4]

        temps = np.unique(dataT)
        alphas = np.zeros_like(temps)
        dalphas = np.zeros_like(temps)
        #dalphas = np.zeros((temps.size,2))
        print(temps)
        # for i,T in enumerate(temps):
        #     alphasT = dataa[dataT==T]
        #     alphas[i] = np.mean(alphasT)
        #     #dalphas[i] = np.std(alphasT, ddof=1)
        #     if len(alphasT>1):
        #         int = interval(alphasT)
        #         alphal,alphau = int[0],int[1]
        #         dalphas[i] = [alphal,alphau]
        #     else:
        #         dalphas[i] = [5e-4, 5e-4]
        # dalphas = dalphas.transpose()
        for i,T in enumerate(temps):
            alphasT = dataa[dataT==T]
            alphas[i] = np.mean(alphasT)
            #dalphas[i] = np.std(alphasT, ddof=1)
            size = alphasT.size
            if size>1:
                ecarttypes = np.std(alphasT, ddof=1)  # écart-types
                size = alphasT.size
                sigma_moyenne = ecarttypes / np.sqrt(size)  # écart-types sur la moyenne
                #############
                ##correction par le facteur de student
                #############
                sigma_corrige = sigma_moyenne * st.t.interval(0.95, size - 1)[1]
                dalphas[i] = sigma_corrige
            else:
                dalphas[i] = dalphas[i-1] if i>0 else 1e-5
            print(size, dalphas[i])
        dalphas = dalphas.transpose()

        print('--------------')
        print(dalphas)

        with plt.style.context(['science', 'scatter']):
            plt.errorbar(temps, alphas*100, xerr=0.05, yerr=dalphas*100, marker='.')


            c = 0.08
            c1 = 5.8/2.4
            T0 = 83.5

            zr = 2.4
            def f_fit(T,c,c1,T0):
                return c*h(1/zr,-zr*c1*(T-T0))

            popt, pcov = curve_fit(f_fit,temps,alphas, p0=[c,c1,T0])
            print(zr,popt)

            TT = np.linspace(75,95,500)

            #plt.plot(TT, [c*h(1/zr,-zr*c1*(T-T0)) for T in TT], 'g--')
            plt.plot(TT, [f_fit(T,*popt)*100 for T in TT], 'g--', label='ajustement pour zr={} cm ($\kappa={}$ %/W, $pdv={}$ K-1)'.format(zr,c,c1))

            zr = 1.4
            def f_fit(T, c, c1, T0):
                return c * h(1 / zr, -zr * c1 * (T - T0))
            popt, pcov = curve_fit(f_fit, temps, alphas, p0=[c, c1, T0])
            print(zr,popt)
            plt.plot(TT, [f_fit(T, *popt)*100 for T in TT], 'k--', label='ajustement pour zr={} cm ($\kappa={}$ %/W, $pdv={}$ K-1)'.format(zr,c,c1))
            plt.legend()
            plt.xlabel("T (°C)")
            plt.ylabel(r"$\alpha$")
            #savefig("alphaT")
        #plt.figure()
    if False:
        plt.scatter(TT,[lpopt(T) for T in TT])
        plt.figure()
        
        plt.plot(bb, [h(0.6,b) for b in bb])
        plt.plot(bb, [h_shift(0.6,b) for b in bb])
        #plt.figure()
    if False:        
        a_arr = np.linspace(0.1,4,40)
        def bopt(a):
            return bb[np.argmax(h(a,bb))]
        fig, ax = plt.subplots(3)
        ax[0].scatter(a_arr,[hmax(a) for a in a_arr])
        ax[0].set_ylabel('h_max')
        ax[1].scatter(a_arr,[bopt(a) for a in a_arr])
        ax[1].set_ylabel('b_opt')
        ax[2].scatter(a_arr,[fwhmb(a)[1]-fwhmb(a)[0] for a in a_arr])
        ax[2].set_ylabel('FWHM en b')
        ax[2].set_xlabel('a')
        plt.figure()
        
        # a_arr2 = np.concatenate((np.linspace(0.1,0.4,20),a_arr))
        # plt.scatter(1/a_arr2,[fwhmb(a)[1]-fwhmb(a)[0] for a in a_arr2])
        # _coeffs = np.polyfit(1/a_arr2,[fwhmb(a)[1]-fwhmb(a)[0] for a in a_arr2],1)
        # print(_coeffs)
        # plt.plot(1/a_arr2, np.polyval(_coeffs,1/a_arr2))
        # plt.figure()
        
        zr_arr = np.linspace(0.1,3,80)    
        plt.plot(zr_arr, [fwhmb(L/2/zr)[1]-fwhmb(L/2/zr)[0] for zr in zr_arr])
        _coeffs = np.polyfit(zr_arr, [fwhmb(L/2/zr)[1]-fwhmb(L/2/zr)[0] for zr in zr_arr],1)
        print(_coeffs)
        plt.plot(zr_arr, np.polyval(_coeffs,zr_arr), label='approximation affine')
        plt.xlabel("$z_R$ (cm)")
        plt.ylabel("$\delta b$ FWHM")
        plt.legend()
        savefig("dbaffine")
        plt.figure()
    if False:
        zr_arr = np.linspace(0.1,3,40)
        markers = ['s','o','^']
        with plt.style.context(['science', 'scatter']):
            for i, lp in enumerate([6.87,6.88,6.9]):
                print(lp)
                plt.plot(zr_arr, [-caract(zr,lp)['dT'] for zr in zr_arr], label='$\Lambda = {}$'.format(lp), marker=markers[i])#, marker = '+')
                plt.legend()
                plt.xlabel("$z_R$ (cm)")
                plt.ylabel("$\delta T$ FWHM (°C)")
                savefig("dt")
    if False:
        plt.plot(TT, [(n2(T)-n1(T))*1e2 for T in TT], label='équation de Sellmeier')
        plt.ylabel(r'$n_2-n_1$ $\left(\times 10^{-2}\right)$')
        plt.xlabel(r'T (°C)')
        #plt.legend()
        #savefig("sellmeier")
        plt.figure()
        zr = 2.4
        bo = bopt(1/zr)
        print(bo)
        plt.plot([2*np.pi / (2*np.pi*(n2(T)-n1(T))/lmbd2 + bo/zr*1e-4) / (1+cexp*(T-T0)) for T in TT], TT)
        #plt.plot([2*np.pi / (2*np.pi*(n2(T)-n1(T))/lmbd2) / (1+cexp*(T-T0)) for T in TT], TT)
        #plt.plot([2*np.pi / (2*np.pi*(n2(T)-n1(T))/lmbd2 + 1.6*1e-4) / (1+cexp*(T-T0)) for T in TT], TT)
        lp_arr = np.linspace(6.75,6.93,100)
        # plt.plot(lp_arr, [caract(zr,lp)['Topt'] for lp in lp_arr]) # pour vérifier validité calculs (not. Tdeb)
        plt.xlabel("$\Lambda$ à 20°C ($\mu$m)")
        plt.ylabel("température optimale (°C)")
        #savefig("topt")
    plt.show()
    
    
