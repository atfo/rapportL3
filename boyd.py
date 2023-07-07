import numpy as np
import scipy.integrate as integr
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from tqdm import tqdm
import sys, argparse


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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', action='store_true')    
    
    zz=np.linspace(0.1,5.1, 41)
    TT = np.linspace(60,120, 201)
    bb = np.linspace(-10,10)
    #plt.plot(bb,[h(200,b) for b in bb])
    #plt.plot(z, [L*alpha(z,100) for z in z])

    zR_arr = np.linspace(0.01, 1.5, 51)
    T_arr = np.linspace(60, 100, 101)

    plt.plot(TT, alpha(1.2, TT,lp))
    plt.figure()

    if parser.parse_args().u:
        im = np.empty((101, 51))
        for i in tqdm(range(101)):
            for j in range(51):
                im[i,j] = alpha(zR_arr[j], T_arr[i],lp)

        print(np.argmax(im))

        np.savez("cmap"+str(lp), im=im)
    else:
        im = np.load("cmap"+str(lp)+".npz")['im']

    plt.imshow(im, extent=(np.amin(zR_arr), np.amax(zR_arr), np.amin(T_arr), np.amax(T_arr)), origin='lower', aspect="auto", cmap=cm.inferno) #, norm=LogNorm())
    plt.xlabel(r'$z_R$ (cm)')
    plt.ylabel('T (°C)')
    plt.colorbar()

    #aa = [alpha(z,T) for z in zz for T in TT]

    #np.savez("vals alpha", alpha=aa, zr=zz, T=TT)

    #plt.imshow(aa)

    #print(opt.minimize(lambda x: -1*alpha(x[0],x[1]),(20,100)))

    #plt.figure()

    plt.show()
