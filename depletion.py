import numpy as np
import matplotlib.pyplot as plt
from odeintw import odeintw
from scipy.integrate import odeint

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

def odeintz(func, z0, t, **kwargs):
    """An odeint-like function for complex valued differential equations."""

    # Disallow Jacobian-related arguments.
    _unsupported_odeint_args = ['Dfun', 'col_deriv', 'ml', 'mu']
    bad_args = [arg for arg in kwargs if arg in _unsupported_odeint_args]
    if len(bad_args) > 0:
        raise ValueError("The odeint argument %r is not supported by "
                         "odeintz." % (bad_args[0],))

    # Make sure z0 is a numpy array of type np.complex128.
    z0 = np.array(z0, dtype=np.complex128, ndmin=1)

    def realfunc(x, t, *args):
        z = x.view(np.complex128)
        dzdt = func(z, t, *args)
        # func might return a python list, so convert its return
        # value to an array with type np.complex128, and then return
        # a np.float64 view of that array.
        return np.asarray(dzdt, dtype=np.complex128).view(np.float64)

    result = odeint(realfunc, z0.view(np.float64), t, **kwargs)

    if kwargs.get('full_output', False):
        z = result[0].view(np.complex128)
        infodict = result[1]
        return z, infodict
    else:
        z = result.view(np.complex128)
        return z

L = 2e-2
coeff = 8.4e-5* 1j # V-1
n1,n2 = 2.1, 2.2

P1 = 0.1
a0 = np.sqrt(6e10*P1)

def func(a,z):
    a1,a2 = a
    return [coeff/n1 * np.conj(a1)*a2, coeff/n2 * (a1**2)]

ai = np.array([a0,0])

zz = np.linspace(-L/2,L/2,100)

a, infodict = odeintz(func,ai,zz, full_output=True)

print(a)
print()
print(np.abs(coeff)*np.abs(a0)**2/(6e10)*L)

plt.plot(zz*1e2, np.abs(a[:,0])**2/(6e7), 'r', label=r'${P}_1$')
plt.plot(zz*1e2, np.abs(a[:,1])**2/(6e7), 'g', label=r"${P}_2$")
plt.plot(zz*1e2, np.square(np.abs(coeff)/n2*np.abs(a0)**2*(zz+L/2))/6e7, 'b--', label=r"${P}_2$ sans déplétion")

plt.xlabel("$z$ (cm)")
plt.ylabel("Puissance (mW)")
plt.legend()
savefig("depl"+str(P1))
plt.show()