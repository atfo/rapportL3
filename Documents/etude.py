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

def idata(fname):
    return np.genfromtxt('./'+fname, delimiter=',', skip_header=1, usemask=True)

def n(T,lmbd):
    a = np.array([5.756, 0.0983, 0.2020, 189.32, 12.52, 1.32e-2])#.reshape(6,1)
    b = np.array([2.860e-6, 4.7e-8, 6.113e-8, 1.516e-4])#.reshape(4,1)
    f = (T-24.5)*(T+570.82)
    p = a[:4] + f*b
    inds = np.sqrt( p[0] + p[1]/(lmbd**2-p[2]**2) + p[3]/(lmbd**2-a[4]**2) - a[5]*lmbd**2)
    return inds

def n1(T):
    return n(T,lmbd1)

def n2(T):
    return n(T,lmbd2)

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
    c= {
        'a' : a,
        'bo' : bo,
        'hm' : hm,
        'bm,bp' : (bm,bp),
        'db' : bp-bm,
        'dkeff_opt' : - bo/zr, # cm-1
        'Topt' : Topt,
        'dT' : dT
        #dT2 = - 2*np.pi * zr /lmbd2 * 1e4 * (n2(100)-n1(100)-n2(80)+n1(80))/20
    }
    return c

caract(2.4,6.9)


if __name__ == "__main__":
    if False:
        depT = idata('dep T basse P.csv')
        dataT = depT[:,3]
        dataa = depT[:,4]
        plt.scatter(dataT, dataa)
        plt.figure()
        
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
    if True:    
        zr_arr = np.linspace(0.1,3,40)
        #with plt.style.context(['science']):#, 'scatter']):
        for lp in [6.85,6.86,6.87,6.88,6.9]:
            print(lp)
            plt.plot(zr_arr, [-caract(zr,lp)['dT'] for zr in zr_arr], label='$\Lambda = {}$'.format(lp))#, marker = '+')
            plt.legend()
            plt.xlabel("$z_R$ (cm)")
            plt.ylabel("$\delta T$ FWHM (°C)")
            savefig("dt")
    plt.show()
    
    
