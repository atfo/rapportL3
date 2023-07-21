#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

import scienceplots
import tikzplotlib
import matplotlib as mpl

from datetime import datetime
import matplotlib.dates as mdates

import pandas as pd

plt.style.use('science')
mpl.rcParams['lines.linewidth'] = 2

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

def idata(fname):
    return np.genfromtxt('./'+fname, delimiter=',', skip_header=1, usemask=True)

if __name__ == "__main__":
#    str2date = lambda x: datetime.strptime(x.decode("utf-8"), '%d/%m/%y %H:%M')
#
#    data = np.genfromtxt('mesc1.csv',dtype=None,names=True, delimiter=',', converters = {0: str2date})
#    print(data.size)
#    print(data[:10])
#    times = data[:][0]
#    print(times)
#    power = data[:][1]
#    #plt.plot(times, power)
#    plt.show()
    
    dateparse = lambda d: datetime.strptime(d, '%d/%m/%y %H:%M')

    #df = pd.read_csv("mesc1.csv", parse_dates=["Heure"], date_parser = dateparse, skip_blank_lines=True)
    df = pd.read_csv("mesc2.csv", skip_blank_lines=True)
    df['Heure'] = pd.to_datetime(df['Heure'], format='%d/%m/%y %H:%M')    
    print(df)
    print(df.info())
    
    times = df["Heure"]
    print(times[:10])

    fig, ax = plt.subplots(1)
    fig.autofmt_xdate()
    plt.scatter(df["Heure"], df["Puissance (W)"], marker='.')

    xfmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    #plt.savefig("mes1c_scatter.pdf", dpi=300)

    fig, ax = plt.subplots(1)
    fig.autofmt_xdate()
    plt.plot(df["Heure"], df["Puissance (W)"])

    xfmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    #plt.savefig("mes1c_plot.pdf", dpi=300)
    savefig("mesc2")
    plt.show()
