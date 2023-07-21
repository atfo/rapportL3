#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy.integrate as integr
import scipy.optimize as opt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from tqdm import tqdm
import sys, argparse
import scienceplots

n = 1 # nbre chemins
g = 0
w = 10
h = 2

# prepare some coordinates
x, y, z = np.indices((100, 30, 30))

# draw cuboids in the top left and bottom right corners, and a link between
# them
cristal = (x < 100) & (y < n*(g+w)+g) & (z < h)

voxelarray = cristal
colors = np.empty(voxelarray.shape + (3,), dtype=object)
blanc = np.ones(3)
colors[cristal] = 0.9*blanc
for i in range(n):
    line = (x < 100) & (y >= i*(g+w)+g) & (y < (i+1)*(g+w)) & (z<8)
    colors[line & (x%(10+i) < (10+i)//2)] = 0.5*blanc
    colors[line & (x%(10+i) >= (10+i)//2)] = 0.7*blanc
    

if __name__ == "__main__":
    # and plot everything
    ax = plt.figure().add_subplot(projection='3d')
    ax.voxels(voxelarray, facecolors=colors, edgecolors=colors)
    ax.set_axis_off()
    plt.savefig("PP.pdf")
    plt.savefig("PP.svg")
    plt.show()

