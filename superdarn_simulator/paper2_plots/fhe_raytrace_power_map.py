#!/usr/bin/env python
import datetime, inspect, os, pickle, shutil, sys
curr_file = inspect.getfile(inspect.currentframe()) # script filename (usually with path)

import glob

import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection

import multiprocessing

from davitpy import pydarn
from davitpy import utils

import mstid
from new_colormaps import viridis

def find_source_file(dt,radar,data_dir='mstid_data/mstid_index'):
    ymd_str         = dt.strftime('%Y%m%d')
    file_pattern    = '{}.*'.format(ymd_str)
    search_str      = os.path.join(data_dir,radar,file_pattern)

    dirs            = glob.glob(search_str)
    for dr in dirs:
        bn          = os.path.basename(dr)
        dt_0        = datetime.datetime.strptime(bn[:13],'%Y%m%d.%H%M')
        dt_1        = datetime.datetime.strptime(bn[14:],'%Y%m%d.%H%M')
        if np.logical_and(dt >= dt_0, dt < dt_1):
            fname   = '{}-{}.p'.format(radar,bn)
            fpath   = os.path.join(dr,fname)
            if os.path.exists(fpath):
                return fpath
    print 'No data found for: {} {!s}'.format(radar,dt)

def pos_hw(x0,y0,x1,y1):
    """
    Convert from x0,y0,x1,y1 format to x0,y0,height,width format.
    """
    return (x0,y0,y1-y0,x1-x0)

def plot_fan(currentData,time,m,cmap_handling='superdarn',cmap=None,plotZeros=True,
        scale=(0,30),param='power',gate_lim=None,beam_lim=None,**kw_args):
    #Figure out which scan we are going to plot...
    if map_time == None:
        timeInx = 0
    else:
        timeInx = (np.where(currentData.time >= map_time))[0]
        if np.size(timeInx) == 0:
            timeInx = -1
        else:
            timeInx = int(np.min(timeInx))

    ngates = np.shape(currentData.data)[2]
    nbeams = np.shape(currentData.data)[1]

    latFull     = currentData.fov.latFull
    lonFull     = currentData.fov.lonFull

    verts = []
    scan  = []
    data  = currentData.data[timeInx,:,:]
    for bm_inx,bm in enumerate(currentData.fov.beams):
        for rg_inx,rg in enumerate(currentData.fov.gates):
#            if goodLatLon[bm,rg] == False: continue

            if np.isnan(data[bm_inx,rg_inx]): continue
            if data[bm_inx,rg_inx] == 0 and not plotZeros: continue
            if gate_lim is not None:
                if gate_lim[0] is not None:
                    if rg  < gate_lim[0]: continue
                if gate_lim[1] is not None:
                    if rg >= gate_lim[1]: continue

            if beam_lim is not None:
                if beam_lim[0] is not None:
                    if bm  < beam_lim[0]: continue
                if beam_lim[1] is not None:
                    if bm >= beam_lim[1]: continue
            scan.append(data[bm_inx,rg_inx])

            x1,y1 = m(lonFull[bm_inx+0,rg_inx+0],latFull[bm_inx+0,rg_inx+0])
            x2,y2 = m(lonFull[bm_inx+1,rg_inx+0],latFull[bm_inx+1,rg_inx+0])
            x3,y3 = m(lonFull[bm_inx+1,rg_inx+1],latFull[bm_inx+1,rg_inx+1])
            x4,y4 = m(lonFull[bm_inx+0,rg_inx+1],latFull[bm_inx+0,rg_inx+1])
            verts.append(((x1,y1),(x2,y2),(x3,y3),(x4,y4),(x1,y1)))

    if cmap_handling == 'matplotlib':
        if cmap is None:
            cmap = matplotlib.cm.jet
        bounds  = np.linspace(scale[0],scale[1],256)
        norm    = matplotlib.colors.BoundaryNorm(bounds,cmap.N)
    elif cmap_handling == 'superdarn':
        colors  = 'lasse'
        cmap,norm,bounds = utils.plotUtils.genCmap(param,scale,colors=colors)

    pcoll = PolyCollection(np.array(verts),edgecolors='face',linewidths=0,closed=False,cmap=cmap,norm=norm,zorder=99)
#    pcoll = PolyCollection(np.array(verts),edgecolors='face',closed=False,cmap=cmap,norm=norm,zorder=99)
    pcoll.set_array(np.array(scan))
    axis    = plt.gca()
    axis.add_collection(pcoll,autolim=False)
    return pcoll

def main(map_time,output_dir='output/map'):
    # Specify which radars we want to include.
    radars = []
    radars.append('fhe')

    # Set figure-wide defaults.
    gate_lim        = (12,50)
    beam_lim        = None
#    data_dir        = 'mstid_data/na_map'
    data_dir        = 'mstid_data/mstid_index'  # MSTID Index is the directory that contains the files
                                                # used for the MSTID classification.

    font            = {'weight':'bold', 'size':10}
    matplotlib.rc('font',**font)

#    map_dict_defaults   = {'lat_0':53.0, 'lon_0':-95, 'width':111e3*60, 'height':111e3*35, 'coords':'geo'}
    map_dict                    =  {'projection':'stere', 'resolution':'l'}
    map_dict['lat_0']           =   44.
    map_dict['lon_0']           =  -93.
    map_dict['width']           =  1550e3
    map_dict['height']          =  1550e3
    map_dict['fillContinents']  = 'None'
    map_dict['fillLakes']       = 'None'
    map_dict_defaults           = map_dict

    fan_dict_defaults   = {}

    # Create a dictionary to hold all of the radar information.
    radar_dct   = {}
    for radar in radars:
        tmp                 = {}
#        source              = find_source_file(map_time,radar,data_dir=data_dir)
#        tmp['source']       = source
        tmp['gate_lim']     = gate_lim
        tmp['beam_lim']     = beam_lim

        radar_dct[radar]    = tmp

    # Set custom beam and gate limits.
#    radar_dct['cve']['beam_lim']        = (0,20)
    radar_dct['fhe']['source']        = 'mstid_data/music_no_mongo/fhe_20121221_paper2_14500kHz_lambda451km/fft/fhe/20121221.1400-20121221.2200/fhe-20121221.1400-20121221.2200.p'


    # Define the subplots
    maps = {}
    tmp = {}
    tmp['radar_dct']                    = radar_dct
    tmp['title']                        = 'Raytrace Simulation\nRaw Ground Scatter Data'
    tmp['dataSet']                      = 'DS000_originalFit'
    tmp['map_dict']                     = map_dict_defaults.copy()
    tmp['fan_dict']                     = fan_dict_defaults.copy()
    tmp['letter']                       = '(a)'
    maps['raw']     = tmp

    tmp = {}
    tmp['radar_dct']                    = radar_dct
    tmp['title']                        = 'Raytrace Simulation\nFiltered Ground Scatter Data'
    tmp['dataSet']                      = 'DS007_detrended'
    tmp['map_dict']                     = map_dict_defaults.copy()
    tmp['fan_dict']                     = fan_dict_defaults.copy()
    tmp['fan_dict']['cmap_handling']    = 'matplotlib'
    tmp['fan_dict']['cmap']             = viridis
    tmp['fan_dict']['scale']            = (-2.,2.)
#    tmp['fan_dict']['cbar_ticks']       = [-2,-1.,0,1.,2.]
    tmp['fan_dict']['cbar_ticks']       = [-1.5,-1.,-0.5,0.,0.5,1.,1.5]
    tmp['letter']                       = '(b)'
    maps['filt']    = tmp

    map_keys        = []
    map_keys.append('raw')
    map_keys.append('filt')

    # Begin plotting loop.
    nr_x            = 1
    nr_y            = len(map_keys)
    fig             = plt.figure(figsize=(nr_x*8,nr_y*6))
    for map_inx,map_key in enumerate(map_keys):
        this_map        = maps[map_key]
        map_dict        = this_map['map_dict']
        fan_dict        = this_map['fan_dict']
        cbar_ticks      = fan_dict.get('cbar_ticks',None)
        letter          = this_map['letter']
        axis            = fig.add_subplot(nr_y,nr_x,map_inx+1)

        axis.set_title(this_map['title'],weight='bold')
        axis.text(0.025,0.9425,map_time.strftime('%d %b %Y %H%M UT'),bbox={'facecolor':'white'},transform=axis.transAxes,size=9,zorder=250)

        m       = utils.plotUtils.mapObj(ax=axis,anchor='SW',dateTime=map_time,showCoords=False,**map_dict)
#        m.drawcoastlines(color='1.',linewidth=4)
        m.drawcoastlines(color='0.755',linewidth=1)
        m.drawmapboundary(fill_color='w')
        for radar in radars:
            source      = this_map['radar_dct'][radar]['source']
            gate_lim    = this_map['radar_dct'][radar].get('gate_lim')
            beam_lim    = this_map['radar_dct'][radar].get('beam_lim')
            beam_lim    = this_map['radar_dct'][radar].get('beam_lim')

            with open(source,'rb') as fl:
                dataObj = pickle.load(fl)

            data_sets   = dataObj.get_data_sets()
            if data_sets == []:
                print 'No data contained in {}!'.format(os.path.basename(source))
                continue

            if this_map['dataSet'] not in data_sets:
                print 'No dataSet {} contained in {}!'.format(this_map['dataSet'],os.path.basename(source))
                continue

            print 'Plotting: {} {} {!s}'.format(os.path.basename(source),this_map['dataSet'],map_time)

            currentData     = getattr(dataObj,this_map['dataSet'])
            rad             = np.array(currentData.metadata['code'],dtype=np.str)[0]
            pcoll           = plot_fan(currentData,map_time,m,beam_lim=beam_lim,gate_lim=gate_lim,**fan_dict)

            pydarn.plotting.overlayRadar(m, rad, zorder=100,fontSize=14)
            pydarn.plotting.overlayFov(m,[rad],dateTime=map_time,
                    rangeLimits=gate_lim,beamLimits=beam_lim,model='GS',
                    zorder=100,lineColor='0.45',lineWidth=1,
                    beams=[7],beamsColors=['gray'])
        m.nightshade(map_time)

        # Figure-wide things.
        cbar_shrink             = 0.85
        cbar_fraction           = 0.150
        cbar_gstext_offset      = -0.125
        cbar_fontsize           = 8
        cbar_gstext_fontsize    = cbar_fontsize

    #    cbar_ax = plt.axes((0.86,0.15,0.020,0.50))
    #    cbar    = fig.colorbar(pcoll,orientation='vertical',cax=cbar_ax)
        cbar    = fig.colorbar(pcoll,orientation='vertical',shrink=cbar_shrink,fraction=cbar_fraction,ax=axis)
        cbar.set_label(r'$\lambda$ Power [dB]',size=9,weight='bold')
        if cbar_ticks is None:
            labels = cbar.ax.get_yticklabels()
            for lbl in labels:
        #        lbl.set_weight('normal')
                lbl.set_size(cbar_fontsize)
            labels[-1].set_visible(False)
        else:
            cbar.set_ticks(cbar_ticks)

        if currentData.metadata.has_key('gscat'):
            if currentData.metadata['gscat'] == 1:
                cbar.ax.text(0.5,cbar_gstext_offset,'Ground\nscat\nonly',ha='center',fontsize=cbar_gstext_fontsize)

        metadata    = currentData.metadata
        txt = 'Coordinates: ' + metadata['coords'] +', Model: ' + metadata['model']
        axis.text(1.001, 1, txt,
                  horizontalalignment='left',
                  verticalalignment='top',
                  rotation='-90',
                  size=8,
                  transform=axis.transAxes)


        axis.text(0.,1.025,letter,transform=axis.transAxes,
                fontdict={'weight':'bold','size':18})

    #axis.set_title(map_time.strftime('%Y %b %m %H%M UT'))
    #fig.tight_layout(pad=2)
    facecolor   = 'w'
    pad_inches  = 0
    filename    = '{}_fhe_raytrace_power_map'.format(map_time.strftime('%Y%m%d_%H%M'))
    filepath    = os.path.join(output_dir,filename)

    extns       = ['png','pdf']
    for extn in extns:
        full_fname  = '{}.{}'.format(filepath,extn)
        fig.savefig(full_fname,facecolor=facecolor,
                bbox_inches='tight',pad_inches=pad_inches)

    plt.close(fig)

        # Add metadata to PDF to know what code generated it.
#            pdf_md = {u'/generating_code':curr_file}
#            mstid.general_lib.pdf_add_metadata(full_fname,pdf_md)
#    if extn == 'pdf':
#        import libxmp
#        import libxmp.utils
#        from libxmp import XMPFiles, consts
#        xmpfile = XMPFiles(file_path=full_fname,open_forupdate=True)
#        xmp     = xmpfile.get_xmp()
#        if xmp is None:
#            xmp = libxmp.XMPMeta()
#
#        xmp.set_property(consts.XMP_NS_DC,u'generating_code',u'hello')
#        res = xmpfile.can_put_xmp(xmp)
#        import ipdb; ipdb.set_trace()
#        xmpfile.close_file()
#
#        from libxmp.utils import file_to_dict,object_to_dict
#        xmp_dict_0  = file_to_dict(full_fname)
#        xmp_dict_1  = object_to_dict(xmp)
#
#        dc = xmp_dict_1[consts.XMP_NS_DC]
#        sys.exit()

if __name__ == '__main__':
    output_dir  = 'output/map'
    mstid.prepare_output_dirs({0:output_dir},clear_output_dirs=False)

    multiproc   = False

#    map_sTime       = datetime.datetime(2014,12,1,18,00)
#    map_eTime       = datetime.datetime(2014,12,1,20,00)

    map_sTime       = datetime.datetime(2012,12,21,16,00)
    map_eTime       = datetime.datetime(2012,12,21,20,00)

    map_times       = []
    curr_time       = map_sTime
    while curr_time < map_eTime:
        map_times.append(curr_time)
        curr_time   = curr_time + datetime.timedelta(minutes=10)

    if multiproc:
        pool = multiprocessing.Pool(6)
        pool.map(main,map_times)
        pool.close()
        pool.join()
    else:
        for map_time in map_times:
            if map_time != datetime.datetime(2012,12,21,16,10):
                continue
            main(map_time)
