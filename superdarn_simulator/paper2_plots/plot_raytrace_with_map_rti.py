#!/usr/bin/env python
import os,sys,shutil
import datetime
import string
import pickle

import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from matplotlib.collections import PolyCollection
import mpl_toolkits.axes_grid.axes_size as Size
from mpl_toolkits.axes_grid import Divider

from davitpy import pydarn

import davitpy
import mstid

import new_colormaps

title_fontdict      = {'weight':'bold',  'size':22}
cbar_label_fontdict = {'weight':'normal','size':22}

plt.rcParams.update({'font.size': 16})
letter_size = 32

def copy_fpath(fpath,copytree=False):
    return # Turn on return when you don't need this function.
    """
    Utility to bring data from original location to this directory for portability.
    """
    print 'Fpath: {}'.format(fpath)
    old_base    = '/data/mstid/statistics/music_scripts/output/raytrace'
    old_path, fname = os.path.split(fpath)

    rel_path    = old_path.lstrip(old_base)

    new_base    = 'raytrace_data'
    new_path    = os.path.join(new_base,rel_path)

    try:
        os.makedirs(new_path)
    except:
        pass
    
    if copytree:
        shutil.rmtree(new_path)
        shutil.copytree(old_path,new_path)
    else:
        try:
            os.remove(os.path.join(new_path,fname))
        except:
            pass
        shutil.copy(fpath,new_path)

    return

def load_edens_objs(radar,panel_time,src_base_path):
    time_str        = panel_time.strftime('%Y%m%d_%H%MUT')
    ret_dict        = {}

    # Load ionos file.
    pkl_fname       = '_'.join([time_str,'250km_edens.p'])
    pkl_fpath       = os.path.join(src_base_path,'edens_alt_p',pkl_fname)
    copy_fpath(pkl_fpath)
    with open(pkl_fpath,'rb') as fl:
        obj = pickle.load(fl)
    ret_dict['edens_alt']       = obj

    return ret_dict

def load_profl_obj(radar,panel_time,beam,src_base_path):
    time_str        = panel_time.strftime('%Y%m%d_%H%MUT')
    # Load generate_ionos_runfile
    
    rdr_str         = '_'.join([time_str,radar,'{:02d}'.format(beam)])
    pkl_fname       = '_'.join([rdr_str,'rto.p'])
    pkl_fpath_0     = os.path.join(src_base_path,'rt_data',rdr_str)
    pkl_fpath       = os.path.join(pkl_fpath_0,pkl_fname)
    copy_fpath(pkl_fpath,copytree=True)
    with open(pkl_fpath,'rb') as fl:
        rto = pickle.load(fl)

    # Required for portability
    rto.outDir = pkl_fpath_0

    print 'Fpath: {}'.format(pkl_fpath)
    return rto

def overlay_fov(radar,map_time,ngates,fov,ax,m,highlight_beam=None):
    if fov is None:
        site    = davitpy.pydarn.radar.radStruct.site(code=radar,dt=map_time)
        fov     = davitpy.pydarn.radar.radFov.fov(site=site,model='IS',coords='geo',ngates=ngates)

    verts = []
    data  = []
    for beam_inx,beam in enumerate(fov.beams):
        if highlight_beam is not None:
            if beam != highlight_beam:
                continue
        for gate_inx,gate in enumerate(fov.gates):
            if beam_inx == len(fov.beams) or \
               gate_inx == len(fov.gates): continue


            rdr_lats = np.zeros((2,2),dtype=np.float32)    
            rdr_lons = np.zeros((2,2),dtype=np.float32)    

            rdr_lats[0,0]   = fov.latFull[beam_inx+0,gate_inx+0]
            rdr_lats[0,1]   = fov.latFull[beam_inx+0,gate_inx+1]
            rdr_lats[1,1]   = fov.latFull[beam_inx+1,gate_inx+1]
            rdr_lats[1,0]   = fov.latFull[beam_inx+1,gate_inx+0]

            rdr_lons[0,0]   = fov.lonFull[beam_inx+0,gate_inx+0]
            rdr_lons[0,1]   = fov.lonFull[beam_inx+0,gate_inx+1]
            rdr_lons[1,1]   = fov.lonFull[beam_inx+1,gate_inx+1]
            rdr_lons[1,0]   = fov.lonFull[beam_inx+1,gate_inx+0]

            if np.any(np.logical_not(np.isfinite([rdr_lats,rdr_lons]))):
                continue

            x1,y1 = m(rdr_lons[0,0], rdr_lats[0,0])
            x2,y2 = m(rdr_lons[0,1], rdr_lats[0,1])
            x3,y3 = m(rdr_lons[1,1], rdr_lats[1,1])
            x4,y4 = m(rdr_lons[1,0], rdr_lats[1,0])

            verts.append(((x1,y1),(x2,y2),(x3,y3),(x4,y4),(x1,y1)))
            
    if highlight_beam is None:
        edgecolors  = '0.1'
        linewidths  = None
    else:
        edgecolors  = 'k'
        linewidths  = 4
    pcoll = PolyCollection(np.array(verts),facecolor=(1,1,1,0),edgecolors=edgecolors,linewidths=linewidths,closed=False,zorder=201)
    ax.add_collection(pcoll,autolim=False)

def plot_map(edens_alt,ax,nel_lim,nel_cmap,radar=None,fov=None,ngates=None,plt_inx=0,
        highlight_beam=None,rasterize_edens=False):


    map_time    = edens_alt['iri_date']
    edens       = np.log10(edens_alt['edens']) # Regular raytrace plots are in log10.
    lons        = edens_alt['lons']
    lats        = edens_alt['lats']
    alt         = edens_alt['alt']

    fig     = ax.get_figure()
    ax_info = {}

    # Wide area North America
    basemap_dict                    =  {'lat_ts':50., 'projection':'stere', 'resolution':'l'}
    basemap_dict['lat_0']           =   60.
    basemap_dict['lon_0']           =  -70.
    basemap_dict['width']           =  6500e3
    basemap_dict['height']          =  6500e3
    basemap_dict['fillContinents']  = 'None'
    basemap_dict['fillLakes']       = 'None'

    m = davitpy.utils.plotUtils.mapObj(ax=ax,datetime=map_time,**basemap_dict)

    # draw parallels and meridians.

    if plt_inx == 0:
#                    [Left,  Right, Top,   Bottom]
        par_labels = [True,  False, False, False]
        mer_labels = [False, False, False, True]
    else:
        par_labels = [False, True,  False, False]
        mer_labels = [False, False, False, True]

    pars = m.drawparallels(np.arange( -80., 81.,20.),color='k',labels=par_labels)
    mers = m.drawmeridians(np.arange(-180.,181.,20.),color='k',labels=mer_labels)

    m.drawcoastlines(color='1.',linewidth=4)
    m.drawcoastlines(linewidth=1)
    m.drawmapboundary(fill_color='w')

    LONS, LATS  = np.meshgrid(lons,lats)
    pcoll       = m.pcolormesh(LONS,LATS,edens,cmap=nel_cmap,vmin=nel_lim[0],vmax=nel_lim[1],
            latlon=True,rasterized=rasterize_edens)

    cbar_label  = r"N$_{el}$ [$\log_{10}(m^{-3})$]"
    ax_info['cbar_pcoll']   = pcoll
    ax_info['cbar_label']   = cbar_label 

    #Overlay radar field of view.
    if radar is not None:
        overlay_fov(radar=radar,map_time=map_time,ngates=ngates,fov=fov,ax=ax,m=m)
        if highlight_beam is not None:
            overlay_fov(radar=radar,map_time=map_time,ngates=ngates,fov=fov,ax=ax,m=m,
                    highlight_beam=highlight_beam)

        davitpy.pydarn.plotting.overlayRadar(m,codes=[radar],markerSize=30,
                fontSize=34,yOffset=-20)

    txt = []
#    txt.append('{0} - Alt: {1:.0f} km'.format(map_time.strftime('%d %b %Y %H%M UT'),float(alt)))
    txt.append('{0} - Alt: {1:.0f} km'.format(map_time.strftime('%d %b %Y %H%M UT'),float(alt)))
    ax.set_title('\n'.join(txt),fontdict=title_fontdict)
    return ax_info

if __name__ == '__main__':
    output_dir = 'output/raytrace'
    mstid.prepare_output_dirs({0:output_dir},clear_output_dirs=True)

    src_base_path   = 'mstid_data/raytrace/fhe_20121221_paper2_14500kHz_lambda451km'

    letter_xpos = -0.10

    nx_ax       = 2
    ny_ax       = 2
    fig_scale   = 2.

    fig = plt.figure(figsize=(fig_scale*6.0,fig_scale*9.))

    nel_lim     = (10,12.25)
#    nel_cmap    = matplotlib.cm.jet
    nel_cmap    = new_colormaps.viridis

    rasterize_edens = True

    ################################################################################
    h_pad       = 0.75
    fig_width   = 0.5

    map_rect    = [0.0, 0.40, 0.8, 0.6]
    rt_rect     = [0.0, 0.20, 0.8, 0.20]
    rti_rect    = [0.0, 0.00, 0.8, 0.15]

    edens_cbar_hgt      = 0.75*(map_rect[3]+rt_rect[3])
    edens_cbar_ypos     = ((map_rect[3]+rt_rect[3])-edens_cbar_hgt)/2. + rt_rect[1]
    edens_cbar_rect     = [0.850, edens_cbar_ypos, 0.175, edens_cbar_hgt]

    rti_cbar_hgt      = 0.90*(rti_rect[3])
    rti_cbar_ypos     = ((rti_rect[3])-rti_cbar_hgt)/2. + rti_rect[1]
    rti_cbar_rect     = [0.850, rti_cbar_ypos, 0.350, rti_cbar_hgt]

    horiz           = [Size.Scaled(1.0)]
    vert            = [Size.Scaled(1.0)]
    map_div         = Divider(fig,map_rect,horiz,vert,aspect=False)
    rt_div          = Divider(fig,rt_rect,horiz,vert,aspect=False)
    rti_div         = Divider(fig,rti_rect,horiz,vert,aspect=False)
    edens_cbar_div  = Divider(fig,edens_cbar_rect,horiz,vert,aspect=False)
    rti_cbar_div    = Divider(fig,rti_cbar_rect,horiz,vert,aspect=False)
    
    pos                 = {}
    pos['a']            = map_div.new_locator( nx=0, ny=0)
    pos['b']            = rt_div.new_locator(  nx=0, ny=0)
    pos['c']            = rti_div.new_locator( nx=0, ny=0)
    pos['edens_cbax']   = edens_cbar_div.new_locator(nx=0, ny=0) 
    pos['rti_cbax']   = rti_cbar_div.new_locator(nx=0, ny=0) 

    ################################################################################
    radar   = 'fhe'
    beam    = 7
    ngates  = 100

    panel_time  = datetime.datetime(2012,12,21,17,00)
    inx         = 0

    # Edensity Map #################################################################
    letter_inx      = inx*ny_ax
    letter          = string.lowercase[letter_inx]

    objs            = load_edens_objs(radar,panel_time,src_base_path)
    edens_alt       = objs['edens_alt']
    
    plot_nr = inx + 1
    ax      = fig.add_subplot(nx_ax,ny_ax,plot_nr)
    ax_info = plot_map(edens_alt,ax,nel_lim=nel_lim,nel_cmap=nel_cmap,
            radar=radar,ngates=ngates,plt_inx=inx,highlight_beam=beam,
            rasterize_edens=rasterize_edens)

    ax.text(letter_xpos,0.95,'({})'.format(letter),fontdict={'size':letter_size,'weight':'bold'},transform=ax.transAxes)
    ax.set_axes_locator(pos[letter])

    pcoll       = ax_info['cbar_pcoll']
    cbar_label  = ax_info['cbar_label']

    cbar        = fig.colorbar(pcoll,orientation='vertical')
    for tl in cbar.ax.get_yticklabels():
        tl.set_size(18)
    cbar.ax.set_axes_locator(pos['edens_cbax'])
    cbar.set_label(cbar_label,fontdict=cbar_label_fontdict,labelpad=15)

    # Plot Edensity Profile ########################################################
    letter_inx      += 1 
    letter          = string.lowercase[letter_inx]

    rto             = load_profl_obj(radar,panel_time,beam,src_base_path)

    rto.readRays()
    rto.readEdens()
    rto.readScatter()

    # Plot rays and electron densities together
    # Plot range markers (every 251 km)
    plot_nr         = nx_ax + inx + 1
    rect            = ''.join([str(x) for x in (ny_ax,nx_ax,plot_nr)])
    ax, aax, cbax   = rto.ionos.plot(panel_time,rect=rect,
                            nel_lim=nel_lim,nel_cmap=nel_cmap,
                            plot_colorbar=False,nel_rasterize=rasterize_edens)
    ax, aax, _      = rto.rays.plot(panel_time, step=10, ax=ax, aax=aax)
    ax.grid()
    time_str        = panel_time.strftime('%d %b %Y %H%M UT')
    txt             = '{} - {} Beam {:.0f}'.format(time_str,radar.upper(),beam)
    ax.set_title(txt,fontdict=title_fontdict)
    ax.text(letter_xpos,0.95,'({})'.format(letter),fontdict={'size':letter_size,'weight':'bold'},transform=ax.transAxes)

    ax.set_axes_locator(pos[letter])

    # Plot RTI Plot ###############################################################
    letter_inx      += 1 
    letter          = string.lowercase[letter_inx]

    #Set plot-wide parameters here.
    yticks          = np.arange(4) * 250. + 250.
    ylim            = (250,1000)

    sDate           = datetime.datetime(2012,12,21,14)
    eDate           = datetime.datetime(2012,12,21,22)
    radar           = 'fhe'
    beam            = 7
    src_file        = 'mstid_data/music_no_mongo/fhe_20121221_paper2_14500kHz_lambda451km/fft/fhe/20121221.1400-20121221.2200/fhe-20121221.1400-20121221.2200.p'
    ds_name         = 'DS000_originalFit'
    cmap            = None
    cmap_handling   = 'superdarn'
    scale           = None
    cbar_ticks      = None
    prefix          = 'Raw Data'

    with open(src_file,'rb') as fl:
        dataObj = pickle.load(fl)

    ax              = fig.add_subplot(111)

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
    res = pydarn.plotting.musicPlot.musicRTI(dataObj,dataSet=ds_name,plotZeros=True,axis=ax, coords='range', yticks=yticks,
            plot_info=False,plot_title=False,ylim=ylim,xBoundaryLimits=xbl,xlim=xlim,secondary_coords='lat',
#            cbar_shrink=cbar_shrink, cbar_fraction=cbar_fraction,cbar_ticks=cbar_ticks,
            cbar_gstext_offset=cbar_gstext_offset, cbar_gstext_fontsize=cbar_gstext_fontsize,beam=beam,plotTerminator=True,
            plot_range_limits_label=False,cmap_handling=cmap_handling,cmap=cmap,scale=scale,plot_cbar=False)

    cbi         = res.cbar_info
    pcoll       = cbi['mappable']
    cbar_label  = cbi['label']
    cbar_ticks  = cbi['ticks']
    
    cbar = fig.colorbar(pcoll,orientation='vertical',aspect=5)
    for tl in cbar.ax.get_yticklabels():
        tl.set_size(18)
    cbar.ax.set_axes_locator(pos['rti_cbax'])
    cbar.set_label(cbar_label,fontdict=cbar_label_fontdict,labelpad=15)
    cbar.set_ticks(cbar_ticks)
    cbar.ax.get_yticklabels()[-1].set_visible(False)
#    cbar.ax.text(0.5,cbar_gstext_offset,'Ground\nscat\nonly',ha='center',fontsize=cbar_gstext_fontsize)

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

    ax.text(letter_xpos,1.15,'({})'.format(letter),fontdict={'size':letter_size,'weight':'bold'},transform=ax.transAxes)
    ax.set_axes_locator(pos[letter])

    ftypes  = []
    ftypes.append('png')
    ftypes.append('pdf')

    for ftype in ftypes:
        fname   = os.path.join(output_dir,'ray_trace_profl.{}'.format(ftype))
        dt_0    = datetime.datetime.now()
        print 'Current time: {!s}'.format(dt_0)
        print '    Saving: {}'.format(fname)
        fig.savefig(fname,bbox_inches='tight')
        dt_1    = datetime.datetime.now()
        dt_t    = dt_1 - dt_0
        print '    Total save time: {!s}'.format(dt_t)

    plt.close('all')
