#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use('science')
import tikzplotlib


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

xx = np.linspace(0,6,1000)

def integr(xx, f):
    S = 0j*np.zeros(xx.size)
    for i in range(1, xx.size):
        S[i] = S[i-1] + (xx[i]-xx[i-1])/2*(f(xx[i]) + f(xx[i-1]))
    return S

def fun(z):
    return 1j*np.exp(-1j*np.pi*z)
def fun2(z):
    return 1j*np.exp(-1j*np.pi*z)*np.sign(np.cos(np.pi*z))
def fun3(z):
    return 1j*np.exp(-1j*np.pi*z)*2/np.pi*np.exp(1j*np.pi*z)
plt.plot(xx,xx)
plt.plot(xx,np.abs(integr(xx,fun)))
plt.plot(xx,np.abs(integr(xx,fun2)))
plt.plot(xx,np.abs(integr(xx,fun3)), linestyle="--")
savefig("QPM")
plt.show()




