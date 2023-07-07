import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

data = np.load("vals alpha.npz")

zz = data['zr']
TT = data['T']
alpha = data['alpha'].reshape(zz.size, TT.size)

print(alpha)

plt.plot(TT, alpha[4,:])

plt.figure()

plt.imshow(alpha, extent=(np.amin(zz), np.amax(zz), np.amin(TT), np.amax(TT)), aspect='auto', cmap=cm.hot)
plt.colorbar()
plt.show()
