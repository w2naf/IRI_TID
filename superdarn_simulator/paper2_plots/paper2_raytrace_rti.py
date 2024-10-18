#!/usr/bin/env python
#This script will create a figure illustrating the processing done in the MUSIC algorithm.
import datetime, inspect, os, pickle, shutil, sys
curr_file = inspect.getfile(inspect.currentframe()) # script filename (usually with path)

import glob
import string

import numpy as np
import scipy as sp
from scipy import stats as stats
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from davitpy import pydarn
from rolling_spect_lib import ConcatenateMusic,MusicFromDataSet

import mstid
from new_colormaps import viridis

import gsmap_geometry

def main(data_dir='mstid_data/mstid_index'):
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 10}

    matplotlib.rc('font', **font)

    #Set plot-wide parameters here.
    output_dir      = 'output/rti'
    mstid.prepare_output_dirs({0:output_dir},clear_output_dirs=True)
    yticks          = np.arange(4) * 250. + 250.
    ylim            = (250,1000)

    data_sources    = []

    tmp = {}
    tmp['radar']        = 'fhe'
    tmp['sDate']        = datetime.datetime(2012,12,21,14)
    tmp['eDate']        = datetime.datetime(2012,12,21,22)
    tmp['beam']         = 7
    tmp['source_files'] = 'mstid_data/music_no_mongo/fhe_20121221_paper2_14500kHz_lambda451km/fft/fhe/20121221.1400-20121221.2200/fhe-20121221.1400-20121221.2200.p'
    tmp['ds_name']      = 'DS000_originalFit'
    tmp['prefix']       = 'Raw Data'
    data_sources.append(tmp)

#    tmp = {}
#    tmp['radar']            = 'fhe'
#    tmp['sDate']            = datetime.datetime(2012,12,21,14)
#    tmp['eDate']            = datetime.datetime(2012,12,21,22)
#    tmp['beam']             = 7
#    tmp['source_files']     = 'mstid_data/music_no_mongo/fhe_20121221_paper2_14500kHz_lambda451km/fft/fhe/20121221.1400-20121221.2200/fhe-20121221.1400-20121221.2200.p'
#    tmp['ds_name']          = 'DS007_detrended'
#    tmp['cmap_handling']    = 'matplotlib'
#    tmp['cmap']             = viridis
#    tmp['scale']            = (-2.0, 2.0)
#    tmp['cbar_ticks']       = np.arange(-1.50,2.00,0.50)
#    tmp['prefix']           = 'Filtered Data'
#    data_sources.append(tmp)

    nx_plots    = 1
    ny_plots    = len(data_sources)
    fig         = plt.figure(figsize=(8*nx_plots,3*ny_plots))

    plot_nr     = 1
    for source in data_sources:
        sDate           = source['sDate']
        eDate           = source['eDate']
        radar           = source['radar']
        beam            = source['beam']
        src_file        = source['source_files']
        ds_name         = source.get('ds_name','DS000_originalFit')
        cmap            = source.get('cmap')
        cmap_handling   = source.get('cmap_handling','superdarn')
        scale           = source.get('scale')
        cbar_ticks      = source.get('cbar_ticks')
        prefix          = source.get('prefix')

        with open(src_file,'rb') as fl:
            dataObj = pickle.load(fl)

        ax              = fig.add_subplot(ny_plots,nx_plots,plot_nr)

        cbar_shrink             = 0.775
        cbar_fraction           = 0.15
        cbar_gstext_offset      = -0.25
        cbar_gstext_fontsize    = 6

        #Plot Original RTI Plot
        currentData = getattr(dataObj,ds_name)
        md          = currentData.metadata
        if md.has_key('gateLimits'): del md['gateLimits']
        if md.has_key('rangeLimits'): del md['rangeLimits']
        if md.has_key('sTime'): del md['sTime']
        if md.has_key('eTime'): del md['eTime']
        xlim        = (sDate, eDate)
        xbl         = (sDate, eDate)
        pydarn.plotting.musicPlot.musicRTI(dataObj,dataSet=ds_name,plotZeros=True,axis=ax, coords='range', yticks=yticks,
                plot_info=False,plot_title=False,ylim=ylim,xBoundaryLimits=xbl,xlim=xlim,secondary_coords='lat',
                cbar_shrink=cbar_shrink, cbar_fraction=cbar_fraction,cbar_ticks=cbar_ticks,
                cbar_gstext_offset=cbar_gstext_offset, cbar_gstext_fontsize=cbar_gstext_fontsize,beam=beam,plotTerminator=False,
                plot_range_limits_label=False,cmap_handling=cmap_handling,cmap=cmap,scale=scale)


        ax.set_ylabel('GS Mapped Geog. Latitude\n'+r'GS Mapped Range $D$ [km]')
        y0,y1   = ax.get_ylim()

        ax_lines    = [2,4,6]
        for ax_line in ax_lines:
            ax_time = sDate + datetime.timedelta(hours=ax_line)
            ax.axvline(ax_time,ls='--',lw=2)

        ax_lines    = [0,2,4,6,8]
        new_xticks  = []
        for ax_line in ax_lines:
            ax_time = sDate + datetime.timedelta(hours=ax_line)
            new_xticks.append(ax_time)
        ax.set_xticks(new_xticks)

        for xtl in ax.get_xticklabels():
            xtl.set_rotation(25)

        ax.set_title(xlim[0].strftime('%d %b %Y'),loc='left')
        title   = []
        title.append(md['name']+' '+'Raytrace Simulation')
        if prefix is not None:
            title.append(prefix)
        ax.set_title('\n'.join(title))
        ax.set_title('Beam %d' % beam,loc='right')

        #Cheat to get lats on detrended data.
        if ds_name == 'DS000_originalFit':
            ytls = []
            for ytl in ax.get_yticklabels():
                ytls.append(ytl.get_text())
        else:
            ax.set_yticklabels(ytls)

        ax.grid()
#        ax.text(-0.15,1.17,'(%s)' % string.lowercase[plot_nr-1],fontdict={'size':16,'weight':'bold'},transform=ax.transAxes)
        plot_nr += 1

    #plt.xticks(rotation=20)
    plt.tight_layout()

    ax  = fig.add_axes([0.9,0.200,0.2,0.70])
    gsmap_geometry.plot_geometry(ax)

    fd      = {'weight':'bold','size':16}
    ypos    = 0.875
    fig.text(0.015,ypos,'(a)',fontdict=fd)
    fig.text(0.900,ypos,'(b)',fontdict=fd)

    file_name   = os.path.join(output_dir,'raytrace_rti')
    fig.savefig(file_name+'.png',bbox_inches='tight')
    fig.savefig(file_name+'.pdf',bbox_inches='tight')
    #fig.clear()

if __name__ == '__main__':
    main()

