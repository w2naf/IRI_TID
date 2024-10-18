#!/usr/bin/env python
import os
import matplotlib
matplotlib.use('Agg')

from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import datetime
from models import *
import utils
import pydarn

import pickle
from struct import unpack

def calculate_scale(data,stddevs=2.,lim='auto'):
    if lim == 'auto':
        mean    = np.mean(np.abs(data))
        std     = np.std(np.abs(data))
        lim     = mean + stddevs*std

    scale   = (-lim,lim)
    ticks   = None
    cmap    = matplotlib.cm.jet

    return scale,ticks,cmap

filename='20141101_1200UT_bks_07.dat'

#with open(filename,'rb') as fl:
#    arr = fl.read()

n_ranges, n_alts = 500, 500

#                self.edens[rtime][raz]['nel'] = array( unpack('{}f'.format(250*250) f.read(250*250*4)) ).reshape((250,250), order='F')
fl = open(filename,'rb')
arr = fl.read(n_ranges*n_alts*4)

res = unpack('{0:d}f'.format(n_ranges*n_alts),arr)
edens = np.array(res).reshape((n_ranges,n_alts),order='F')


range_step  = 5
alt_step    = 1
ranges  = np.arange(n_ranges)*5
alts    = np.arange(60,560)
################################################################################
fig  = plt.figure(figsize=(10,6))
ax   = fig.add_subplot(111)

verts = []
data  = []
for range_inx,rnge in enumerate(ranges):
    for alt_inx,alt in enumerate(alts):
        x1,y1 = rnge-range_step/2.,alt-alt_step/2.
        x2,y2 = rnge+range_step/2.,alt-alt_step/2.
        x3,y3 = rnge+range_step/2.,alt+alt_step/2.
        x4,y4 = rnge-range_step/2.,alt+alt_step/2.
        verts.append(((x1,y1),(x2,y2),(x3,y3),(x4,y4),(x1,y1)))
        data.append(edens[range_inx,alt_inx])

data = np.array(data)
scale,ticks,cmap = calculate_scale(data)
bounds  = np.linspace(scale[0],scale[1],256)
norm    = matplotlib.colors.BoundaryNorm(bounds,cmap.N)

pcoll = PolyCollection(np.array(verts),edgecolors='face',linewidths=0,closed=False,cmap=cmap,norm=norm)
pcoll.set_array(data)
ax.add_collection(pcoll,autolim=False)

ax.set_xlim(ranges.min(), ranges.max())
ax.set_ylim(alts.min(), alts.max())

ax.set_xlabel('Range [km]')
ax.set_ylabel('Altitude [km]')

cbar    = fig.colorbar(pcoll,orientation='vertical',shrink=0.60,pad=.10,ticks=ticks)
txt     = r'IRI Electron Density [m$^{-3}$]'
cbar.set_label(txt)

txt = []
txt.append('IRI Electron Density')
txt.append(filename)
ax.set_title('\n'.join(txt))

fig.tight_layout()
fig.savefig(os.path.join('data','{0}.png'.format(filename)),bbox_inches='tight')
