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
import scienceplots
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

L=2 #cm
lp = 6.874 #um
lmbd1 = 1.064 #um
lmbd2 = 0.532

#ref: gayer2008
def n(T,lmbd):
    a = np.array([5.756, 0.0983, 0.2020, 189.32, 12.52, 1.32e-2]).reshape(6,1)
    b = np.array([2.860e-6, 4.7e-8, 6.113e-8, 1.516e-4]).reshape(4,1)
    f = (T-24.5)*(T+570.82)
    p = a[:4] + f*b
    inds = np.sqrt( p[0] + p[1]/(lmbd**2-p[2]**2) + p[3]/(lmbd**2-a[4]**2) - a[5]*lmbd**2)
    return inds

def n1(T):
    return n(T,lmbd1)

def n2(T):
    return n(T,lmbd2)

def complex_quad(func, a, b, **kwargs):
    def real_func(x):
        return np.real(func(x))
    def imag_func(x):
        return np.imag(func(x))
    real_integral = integr.quad(real_func, a, b, **kwargs)
    imag_integral = integr.quad(imag_func, a, b, **kwargs)
    return real_integral[0] + 1j*imag_integral[0]

def h(a,b):
    def f(t):
        return np.exp(1j*b*t)/(1+1j*t)
    #return (np.abs(complex_quad(f,-a,a, limit=1000)))**2 / 4*a
    N = 1000
    x_array = np.linspace(-a, a, N).reshape(N,1)
    f_array = np.array([f(x) for x in x_array])
    f_array = f_array.reshape(N, f_array.size // N)
    #print(f_array.shape)
    integral = np.trapz(f_array, x_array, axis=0)
    return 1 / (4*a) * np.square(np.abs(integral))

def alpha(zr,T,lp):
    n1=n(T,lmbd1)
    n2=n(T,lmbd2)
    a=L/(2*zr)
    #print('a: {}'.format(a))
    dk=2*np.pi*(2*n1/lmbd1-n2/lmbd2+1/lp)*1e4 #de um-1 a cm-1
    b = dk*zr
    #print(dk)
    #print('b: {}'.format(b))
    return L*h(a,b)/n1/n2

#print(h(1,2))
#print(alpha(10,100))

def deltak(T):
    n1=n(T,lmbd1)
    n2=n(T,lmbd2)
    return 2*np.pi*(2*n1/lmbd1-n2/lmbd2)*1e4 #de um-1 a cm-1

def grad(f,x):
    return (np.roll(f,1)-f)/(np.roll(x,1) - x)

def savefig(name):
    tpl_fix(plt.gcf())
    tikzplotlib.save(name+".tex")

if __name__ == "__main__":
    TT = np.linspace(40,110,500)
    plt.plot(TT, deltak(TT))
    plt.ylabel(r'$\Delta k$ ($cm^{-1}$)')
    plt.xlabel(r'T (°C)')
    plt.figure()
    def dT(TT):
        beta = 15e-6
        n2a = n2(TT)
        n1a = n1(TT)
        #grads = np.gradient(n2,TT)-np.gradient(n1,TT)
        grads = grad(n2a,TT) - grad(n1a,TT) 
        return 0.4429*lmbd1/L*1e-4 / (grads - beta*(n2a-n1a))
    plt.plot(TT, dT(TT), '.')
    plt.ylabel(r'$\delta T$ (°C)')
    plt.xlabel(r'T (°C)')
    plt.figure()
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', action='store_true')    
    s1 = 251
    s2 = 201
    a_arr = np.linspace(0.01, 8, s1)
    b_arr = np.linspace(-4, 4, s2)
    
    bb = np.linspace(-4, 4, 301)
    for a in [0.5,1,2,2.886,3,100]:
        plt.plot(bb, h(a, bb), label='a={}'.format(a))
        #plt.plot(bb, a*(np.sinc(a*bb/np.pi))**2, label='approx onde plane')
    plt.xlabel(r'$b$')
    plt.ylabel(r'$h(a,b)$')
    plt.legend()
    #plt.savefig("bk-factor.pdf", dpi=300)
    tpl_fix(plt.gcf())
    tikzplotlib.save("bk-factor.tex")
    plt.figure()

    if parser.parse_args().u:
        im = np.empty((s2, s1))
        for i in tqdm(range(s2)):
            for j in range(s1):
                im[i,j] = h(a_arr[j], b_arr[i])
        np.savez("h_cmap"+str(lp), im=im)
    else:
        im = np.load("h_cmap"+str(lp)+".npz")['im']
    

    ind = (np.unravel_index(im.argmax(), im.shape))
    print(b_arr[ind[0]],a_arr[ind[1]])
    plt.imshow(im, extent=(np.amin(a_arr), np.amax(a_arr), np.amin(b_arr), np.amax(b_arr)), origin='lower', aspect="auto", cmap=cm.inferno) #, norm=LogNorm())
    plt.xlabel(r'$a$')
    plt.ylabel(r'$b$')
    plt.colorbar()
    #plt.savefig("bk-factor-cmap.pdf", dpi=300)    
    savefig("bk-factor-cmap")    
    plt.show()

