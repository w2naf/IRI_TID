#!/usr/bin/env python
import os,sys,shutil
import datetime

import pickle
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import string

import numpy as np
import davitpy
import mstid

output_dir = 'output/raytrace'
mstid.prepare_output_dirs({0:output_dir},clear_output_dirs=True)

p_files = []
p_files.append('/data/mstid/statistics/music_scripts/output/raytrace/fhe_20121221_paper2_14500kHz_lambda451km/rt_data/20121221_1700UT_fhe_07/20121221_1700UT_fhe_07_rtrun.p')
p_files.append('/data/mstid/statistics/music_scripts/output/raytrace/fhe_20121221_paper2_14500kHz_lambda451km/rt_data/20121221_1720UT_fhe_07/20121221_1720UT_fhe_07_rtrun.p')

nx_ax   = 1
ny_ax   = len(p_files)

fig = plt.figure(figsize=(10,6.5))
plt.rcParams.update({'font.size': 14})

for inx,p_file in enumerate(p_files):

    with open(p_file,'rb') as fl:
        rto = pickle.load(fl)

    rto.readRays()
    rto.readEdens()
    rto.readScatter()

    rt_time = rto.time[0]
    radar   = rto.radar.code[0]
    beam    = rto.beam

    # Plot rays and electron densities together
    # Plot range markers (every 250 km)
    rect    = ''.join([str(x) for x in (ny_ax,nx_ax,inx+1)])
    ax, aax, cbax = rto.ionos.plot(rt_time,rect=rect)
    ax, aax, cbax = rto.rays.plot(rt_time, step=10, ax=ax, aax=aax)
    ax.grid()
    time_str    = rt_time.strftime('%d %b %Y %H%M UT')
    txt         = 'Raytrace: {} Beam {:.0f} [{}]'.format(radar.upper(),beam,time_str)
    ax.set_title(txt)
    ax.text(-0,1.00,'(%s)' % string.lowercase[inx],fontdict={'size':16,'weight':'bold'},transform=ax.transAxes)

fig.tight_layout()

ftypes  = []
ftypes.append('png')
ftypes.append('pdf')

for ftype in ftypes:
    fname = os.path.join(output_dir,'ray_trace_profl.{}'.format(ftype))
    fig.savefig(fname,bbox_inches='tight')

plt.close('all')
