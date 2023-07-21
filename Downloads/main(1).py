import numpy as np
import scipy.integrate as int
import scipy.optimize as opt
import matplotlib.pyplot as plt

a = np.array([5.756, 0.0983, 0.2020, 189.32, 12.52, 1.32e-2])
b = np.array([2.860e-6, 4.7e-8, 6.113e-8, 1.516e-4])

L=2 #cm
lp = 6.83 #um
lmbd1 = 1.064 #um
lmbd2 = 0.532

def n(T,lmbd):
    f = (T-24.5)*(T+570.82)
    p = a[:4] + f*b
    return np.sqrt( p[0] + p[1]/(lmbd**2-p[2]**2) + p[3]/(lmbd**2-a[4]**2) - a[5]*lmbd**2)

def complex_quad(func, a, b, **kwargs):
    def real_func(x):
        return np.real(func(x))
    def imag_func(x):
        return np.imag(func(x))
    real_integral = int.quad(real_func, a, b, **kwargs)
    imag_integral = int.quad(imag_func, a, b, **kwargs)
    return real_integral[0] + 1j*imag_integral[0]

def h(a,b):
    def f(t):
        return np.exp(-1j*b*t)/(1+1j*t)
    return (np.abs(complex_quad(f,-a,a, limit=500)))**2 / 4*a

def alpha(zr,T):
    n1=n(T,lmbd1)
    n2=n(T,lmbd2)
    a=L/(2*zr)
    print('a: {}'.format(a))
    dk=2*np.pi*(2*n1/lmbd1-n2/lmbd2-1/lp)*1e4 #de um-1 a cm-1
    b = dk*zr
    print(dk)
    print('b: {}'.format(b))
    return h(a,b)/n1/n2

print(h(1,2))
print(alpha(10,100))

z=np.linspace(-100,100)
bb = np.linspace(-10,10)
#plt.plot(bb,[h(200,b) for b in bb])
plt.plot(z, [L*alpha(z,100) for z in z])

print(opt.minimize(lambda x: -1*alpha(x[0],x[1]),(20,100)))

#plt.figure()

plt.show()
