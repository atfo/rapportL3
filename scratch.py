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

plt.style.use('science')

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
    return 0.143*L*h(a,b)/n1/n2

#print(h(1,2))
#print(alpha(10,100))

def deltak(T):
    n1=n(T,lmbd1)
    n2=n(T,lmbd2)
    return 2*np.pi*(2*n1/lmbd1-n2/lmbd2)*1e4 #de um-1 a cm-1

def grad(f,x):
    return (np.roll(f,1)-f)/(np.roll(x,1) - x)

Lc = np.array([6.83,6.86,6.9,6.93,6.96])
Tc = np.array([103,88,67,50,33])

def edind(Lc):
    return 0.532*(1/Lc)# - 1.6e-4/(2*np.pi))

def dn(T):
    coeffs = np.polyfit(Tc-80, edind(Lc), deg=2)
    return np.polyval(coeffs, T-80)

def lpopt(T):
    return 2*np.pi / (1.6e-4 + 2*np.pi*dn(T)/lmbd2)

def alphaC(zr,T,lp):
    a = L/(2*zr)
    b = zr * 2*np.pi* (1/lp - dn(T)/lmbd2 ) * 1e4
#    print('a={}, b={}'.format(a,b))
    return 0.143*L/n1(T)/n2(T) * h(a,b)

cexp = 1.5e-5

def alpha_corr(zr,T,lp):    
    LT = L*(1+cexp*(T-20))
    lpT = lp*(1+cexp*(T-20))
    a = LT/(2*zr)
    b = zr * 2*np.pi* (1/lpT - dn(T)/lmbd2 ) * 1e4
    return 0.143*LT/n1(T)/n2(T) * h(a,b)

if __name__ == "__main__":
    TT = np.linspace(30,150,500)
    plt.plot(TT, deltak(TT))
    plt.ylabel(r'$\Delta k$ ($cm^{-1}$)')
    plt.xlabel(r'T (°C)')
    plt.figure()
    # n2-n1 from Covesion data  
    plt.plot(Tc, edind(Lc), '+', label='données Covesion')
    plt.plot(TT, n2(TT)-n1(TT), '--', label='équation de Sellmeier')
    coeffs = np.polyfit(Tc-80, edind(Lc), deg=2)
    print(coeffs)
    plt.plot(TT, np.polyval(coeffs,TT-80), '--', label='ajustement quadratique')
    plt.xlabel(r'T (°C)')
    plt.ylabel(r'$n_2-n_1$')
    plt.xlabel(r'T (°C)')
    plt.legend()
    plt.figure()
    plt.plot(TT, lpopt(TT), label='sans dilatation thermique')
    plt.plot(TT, lpopt(TT)*(1-cexp*(TT-20)), label='avec dilatation thermique')
    plt.ylabel(r"période d'inversion $\Lambda$ optimale ($\mu$m)")
    plt.xlabel(r'T (°C)')
    plt.legend()
    plt.figure()
    # h(T) from C data
    plt.plot(TT, alphaC(0.35,TT,6.9), label='sans dilatation thermique')
    plt.plot(TT, [alpha_corr(0.35,T,6.9) for T in TT], label='avec dilatation thermique')
    plt.ylabel(r'$\alpha$(T) (W/W$^2$)')
    plt.xlabel(r'T (°C)')
    plt.legend()
    plt.figure()

    # temperature bandwidth

    def dT(TT):
        beta = 15e-6
        n2a = n2(TT)
        n1a = n1(TT)
        #grads = np.gradient(n2,TT)-np.gradient(n1,TT)
        grads = grad(n2a,TT) - grad(n1a,TT) 
        return 0.4429*lmbd1/L*1e-4 / (grads - beta*(n2a-n1a))
#    plt.plot(TT, dT(TT), '.')
#    plt.ylabel(r'$\delta T$ (°C)')
#    plt.xlabel(r'T (°C)')
#    plt.figure()


    parser = argparse.ArgumentParser()
    parser.add_argument('-u', action='store_true')    
    s1 = 201
    s2 = 101
    a_arr = np.linspace(0.01, 8, s1)
    b_arr = np.linspace(-4, 4, s2)

    if parser.parse_args().u:
        im = np.empty((s2, s1))
        for i in tqdm(range(s2)):
            for j in range(s1):
                im[i,j] = h(a_arr[j], b_arr[i])
        np.savez("scratch_cmap"+str(lp), im=im)
    else:
        im = np.load("scratch_cmap"+str(lp)+".npz")['im']
    

    ind = (np.unravel_index(im.argmax(), im.shape))
    print(b_arr[ind[0]],a_arr[ind[1]])
    plt.imshow(im, extent=(np.amin(a_arr), np.amax(a_arr), np.amin(b_arr), np.amax(b_arr)), origin='lower', aspect="auto", cmap=cm.inferno) #, norm=LogNorm())
    plt.xlabel(r'$a$')
    plt.ylabel(r'$b$')
    plt.colorbar()
    plt.show()

