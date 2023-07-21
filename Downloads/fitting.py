#!/usr/bin/python3
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


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


# In[6]:


L=2 #cm
#lp = 6.874 #um
lmbd1 = 1.064 #um
lmbd2 = 0.532

print("n1={} et n2={} à 100°C".format(n1(10), n2(100)))
TT = np.linspace(30,150,500)

# In[8]:


cexp = 1.5e-5

def alpha(zr,T,lp,c):  
    LT = L*(1+cexp*(T-20))
    lpT = lp*(1+cexp*(T-20))
    a = LT/(2*zr)
    b = zr * 2*np.pi* (1/lpT - (n2(T)-n1(T))/lmbd2 ) * 1e4
    return c*LT/n1(T)/n2(T) * h(a,b)


# In[9]:


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


# In[14]:


bb=np.linspace(-20,20,1000)
a=0.4
plt.plot(bb,h(a,bb),label='a={}'.format(a))
plt.legend()
m = np.max(h(a,bb))
print(m)
plt.hlines(m,-20,20)
plt.hlines(m/2,-20,20)
plt.show()


# In[ ]:




