import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.constants import physical_constants
from tqdm import tqdm
from scipy.optimize import curve_fit

from tikzplotlib import save as tikz_save, clean_figure as tikz_clean

eta = np.sqrt(physical_constants["vacuum mag. permeability"][0] / physical_constants["vacuum electric permittivity"][0])
L = 40e-3
lambda1 = 1550.1198e-9
lambda2 = 1050.3399e-9
lambda3 = 626.082e-9
lambda3 = 1/(1/lambda1 + 1/lambda2)
d33 = 25e-12

def delta_k(n_array):
    return 2 * np.pi * (n_array[2]/lambda3 - n_array[1]/lambda2 - n_array[0]/lambda1)

def h(a,b,c):
    N = 1000
    f = lambda x: np.exp(-1j*b*x) / ((1 + 1j*x) * (1 + 1j*c*x))
    x_array = np.linspace(-a, a, N)
    f_array = np.array([f(x) for x in x_array])
    integral = np.trapz(f_array, x_array)
    
    return 1 / (4*a) * np.square(np.abs(integral))

def a(zR):
    return L / (2*zR)

def b(zR, n_array, Lambda):
    return (delta_k(n_array) - 2*np.pi/Lambda) * zR

def c(zR, n_array, Lambda):
    w_eff_squared = 1 / (np.pi / zR * (n_array[0]/lambda1 + n_array[1]/lambda2 + n_array[2]/lambda3))
    return (b(zR, n_array, Lambda) + np.pi*L/(a(zR)*Lambda)) * 1/(np.pi*L) * 1/(n_array[0]/lambda1 + n_array[1]/lambda2 + n_array[2]/lambda3)

def Z0(n_array):
    return 32 * np.pi**2 * eta / (lambda1 * lambda2 * (n_array[0]/lambda1 + n_array[1]/lambda2 + n_array[2]/lambda3)**2)

def n(T, lambda_i): # T in CELSIUS !!!
    f = lambda T: (T - 25.4) * (T + 570.82)
    a1 = 5.756
    a2 = 0.0983
    a3 = 0.2020
    a4 = 189.32
    a5 = 12.52
    a6 = 1.32e-2
    b1 = 2.86e-6
    b2 = 4.7e-8
    b3 = 6.113e-8
    b4 = 1.516e-4
    
    return np.sqrt(a1 + b1*f(T) + (a2 + b2*f(T))/(lambda_i**2/(1e-6)**2 - (a3 + b3*f(T))**2) + (a4 + b4*f(T))/(lambda_i**2/(1e-6)**2 - a5**2) - a6*lambda_i**2/(1e-6)**2)

def alpha(T, zR, Lambda):
    n_array = [n(T, lambda1), n(T, lambda2), n(T, lambda3)]
    return Z0(n_array) * L * (2*d33/np.pi)**2 / lambda3**3 * h(a(zR), b(zR, n_array, Lambda), c(zR, n_array, Lambda))

#%% 2D plot
zR_arr = np.linspace(2, 12, 100)*1e-3
T_arr = np.linspace(190, 230, 100)

im = np.empty((100, 100))
for i in tqdm(range(100)):
    for j in range(100):
        im[i,j] = alpha(T_arr[i], zR_arr[j], 11.12e-6)
        
plt.imshow(im, extent=(2, 12, 230, 190), aspect="auto", cmap=cm.hot)

plt.imsave("alpha_T_zR_2D.svg", im)

plt.colorbar()

plt.show()


#%% 1D plot

zR_arr = np.array([3, 7, 11])*1e-3
T_arr = np.linspace(170, 210, 500)
alpha_arr = []
for zR in zR_arr:
    alpha_arr.append([alpha(T, zR, 11.12e-6) for T in T_arr])

plt.figure()
for i in range(len(zR_arr)):
    plt.plot(T_arr, alpha_arr[i], label=rf"$z_R = ${zR_arr[i]*1e-3} mm")
plt.xlabel("Température du cristal (\si{\celsius})")
plt.ylabel(r"$\alpha_{th}$")
plt.legend()

# tikz_clean()
# tikz_save("alpha_T_zR_1D.tikz")

#%% delta_k in case of unpolarised non-linear medium

T_array = np.linspace(0, 200, 500)
delta_k_array = np.array([delta_k([n(T, lambda1), n(T, lambda2), n(T, lambda3)]) for T in T_array])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(T_array, delta_k_array)
ax.set_ylabel("$\Delta k$ (\si{\per\meter})")
ax2 = ax.twinx()
ax2.plot(T_array, 2*np.pi / delta_k_array * 1e6, c="C1")
ax2.set_xlabel("Température du cristal (\si{\celsius})")
ax2.set_ylabel("$2\pi / \Delta k$ (\si{\micro\meter})")

# tikz_save("delta_k.tikz")


# plt.plot(T_array, delta_k_array)
# plt.ylabel("$\Delta k$")
# plt.twinx()
# plt.plot(T_array, 2*np.pi / delta_k_array)
# plt.xlabel("Température du cristal (\si{\celsius})")
# plt.ylabel("$\Delta k$")
# tikz_clean()
# tikz_save("collimated_h.tikz")
    
#%% plots h(a,b,c)


x_array = np.linspace(-10, 10, 500)
y_array = [h(1, x, 0) for x in x_array]
    

plt.figure()
plt.subplots_adjust(wspace=0.5)
plt.subplot(121)
plt.plot(x_array, y_array)
plt.xlabel("$b$")
plt.ylabel("$h(1, b, 0)$")

tikz_clean()
tikz_save("focalised_h.tikz")

y_array = [h(200, x, 0) for x in x_array]

plt.subplot(122)
plt.plot(x_array, y_array)
plt.xlabel("$b$")
plt.ylabel("$h(200, b, 0)$")
# tikz_clean()
# tikz_save("collimated_h.tikz")


#%% Experimental data fitting

T_array = np.array([170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 187.5, 188, 188.5, 189, 189.25, 189.5, 189.75, 190, 190.25, 190.5, 190.75, 191, 191.5, 192, 192.5, 193, 194, 195, 196, 197, 198, 199, 200])
T_array = T_array[4:]

P_array = np.array([0.070, 0.071, 0.071, 0.072, 0.071, 0.076, 0.074, 0.077,
                    0.075, 0.077, 0.079, 0.075, 0.085, 0.082, 0.090, 0.097,
                    0.098, 0.149, 0.330, 0.206, 0.725, 3.52, 7.4, 15.3, 24.2,
                    28.8, 23.2, 12.3, 3.06, 0.311, 3.73, 0.602, 0.843, 0.332,
                    0.218, 0.191, 0.161, 0.137, 0.1028, 0.134, 0.106])
P_array = P_array[4:]


# T_array = np.array([180, 181, 182, 183, 184, 185, 186, 187, 187.5, 188, 188.5, 189, 189.25, 189.5, 189.75, 190, 190.25, 190.5, 190.75, 191, 191.5, 192, 192.5, 193, 194, 195, 196, 197, 198, 199, 200], dtype=np.float64)

# P_array = np.array([0.079, 0.075, 0.085, 0.082, 0.090, 0.097,
#                     0.098, 0.149, 0.330, 0.206, 0.725, 3.52, 7.4, 15.3, 24.2,
#                     28.8, 23.2, 12.3, 3.06, 0.311, 3.73, 0.602, 0.843, 0.332,
#                     0.218, 0.191, 0.161, 0.137, 0.1028, 0.134, 0.106], dtype=np.float64)


def f(T, zR, Lambda, amp):
    n_array = [n(T, lambda1), n(T, lambda2), n(T, lambda3)]
    _a = a(zR)
    _b = b(zR, n_array, Lambda)
    _c = c(zR, n_array, Lambda)
    
    return amp * h(_a, _b, _c)

def f_fit(T, zR, Lambda, amp):
    res = []
    for t in T:
        n_array = [n(t, lambda1), n(t, lambda2), n(t, lambda3)]
        _a = a(zR*1e-3)
        _b = b(zR*1e-3, n_array, Lambda*1e-6)
        _c = c(zR*1e-3, n_array, Lambda*1e-6)
        
        res.append(amp * h(_a, _b, _c))
    return np.array(res)

popt, pcov = curve_fit(f_fit, T_array, P_array, p0=[7.1, 11.145, 26], bounds=([0, 0, 0],[np.inf, np.inf, np.inf]))
print(popt)

#%% Experimental data plotting

T_array_plot = np.linspace(T_array[0], T_array[-1], 200)

f_array = [f(T, popt[0]*1e-3, popt[1]*1e-6, popt[2]) for T in T_array_plot]
f_th_array = [f(T, 7*1e-3, 11.17*1e-6, popt[2]) for T in T_array_plot]

plt.figure()
plt.scatter(T_array, P_array, marker="o")
plt.plot(T_array_plot, f_array, color="black")
plt.plot(T_array_plot, f_th_array, color="black", ls="--")

plt.xlabel(r"$T$ (\si{\celsius})")
plt.ylabel(r"Power (\si{\milli\watt})")

plt.tick_params(direction="in")
# tikz_clean()
# tikz_save("puissance_temp.tikz")

