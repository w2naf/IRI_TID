#!/usr/bin/env python3
"""
Figure_2b_TID_RayTrace.py
Nathaniel A. Frissell
February 2024

This script is used to generate Figure 2b of the Frissell et al. (2024)
GRL manuscript on multi-instrument measurements of AGWs, MSTIDs, and LSTIDs.
"""
import os 

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib
mpl = matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import astropy.units as ap_u
from astropy.time import Time as ap_Time
from astropy.coordinates import get_sun as ap_get_sun
from astropy.coordinates import EarthLocation as ap_EarthLocation
from astropy.coordinates import AltAz as ap_AltAz

from matplotlib.transforms import Affine2D, Transform
import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.projections import polar
from mpl_toolkits.axisartist.grid_finder import FixedLocator, DictFormatter

from pylap.raytrace_2d import raytrace_2d 

mpl.rcParams['font.size']               = 16.0
mpl.rcParams['axes.labelsize']          = 'xx-large'
mpl.rcParams['axes.titlesize']          = 'xx-large'
mpl.rcParams['figure.titlesize']        = 'xx-large'
mpl.rcParams['legend.fontsize']         = 'large'
mpl.rcParams['legend.title_fontsize']   = None

mpl.rcParams['xtick.labelsize']         = 32
mpl.rcParams['ytick.labelsize']         = 32

title_fontsize      = 28

def curvedEarthAxes(rect=111, fig=None, minground=0., maxground=2000, minalt=0,
                    maxalt=1500, Re=6371., scale_heights=1.,nyticks=4, nxticks=4):
    """Create curved axes in ground-range and altitude

    Parameters
    ----------
    rect : Optional[int]
        subplot spcification
    fig : Optional[pylab.figure object]
        (default to gcf)
    minground : Optional[float]

    maxground : Optional[int]
        maximum ground range [km]
    minalt : Optional[int]
        lowest altitude limit [km]
    maxalt : Optional[int]
        highest altitude limit [km]
    Re : Optional[float] 
        Earth radius in kilometers
    nyticks : Optional[int]
        Number of y axis tick marks; default is 5
    nxticks : Optional[int]
        Number of x axis tick marks; deafult is 4

    Returns
    -------
    ax : matplotlib.axes object
        containing formatting
    aax : matplotlib.axes object
        containing data

    Example
    -------
        import numpy as np
        ax, aax = curvedEarthAxes()
        th = np.linspace(0, ax.maxground/ax.Re, 50)
        r = np.linspace(ax.Re+ax.minalt, ax.Re+ax.maxalt, 20)
        Z = exp( -(r - 300 - ax.Re)**2 / 100**2 ) * np.cos(th[:, np.newaxis]/th.max()*4*np.pi)
        x, y = np.meshgrid(th, r)
        im = aax.pcolormesh(x, y, Z.T)
        ax.grid()

    Nathaniel A. Frissell, February 2024
    Adapted from code by Sebastien de Larquier, April 2013
    """
    ang         = maxground / Re
    minang      = minground / Re
    angran      = ang - minang
    angle_ticks = [(0, "{:.0f}".format(minground))]
    while angle_ticks[-1][0] < angran:
        tang = angle_ticks[-1][0] + 1./nxticks*angran
        angle_ticks.append((tang, "{:.0f}".format((tang-minang)*Re)))

    grid_locator1   = FixedLocator([v for v, s in angle_ticks])
    tick_formatter1 = DictFormatter(dict(angle_ticks))

    altran      = float(maxalt - minalt)
    alt_ticks   = [(minalt+Re, "{:.0f}".format(minalt))]
    while alt_ticks[-1][0] <= Re+maxalt:
        crd_alt  = altran / float(nyticks) + alt_ticks[-1][0]
        real_alt = (crd_alt - Re) / scale_heights

        alt_ticks.append((crd_alt, "{:.0f}".format(real_alt)))
    _ = alt_ticks.pop()
    grid_locator2   = FixedLocator([v for v, s in alt_ticks])
    tick_formatter2 = DictFormatter(dict(alt_ticks))

    tr_rotate       = Affine2D().rotate(np.pi/2-ang/2)
    tr_shift        = Affine2D().translate(0, Re)
#    tr              = polar.PolarTransform(apply_theta_transforms=False) + tr_rotate
    tr              = polar.PolarTransform() + tr_rotate

    grid_helper = \
        floating_axes.GridHelperCurveLinear(tr, extremes=(0, angran, Re+minalt,
                                                          Re+maxalt),
                                            grid_locator1=grid_locator1,
                                            grid_locator2=grid_locator2,
                                            tick_formatter1=tick_formatter1,
                                            tick_formatter2=tick_formatter2,)

    if not fig: fig = plt.gcf()
    ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)

    # adjust axis
    ax1.axis["left"].label.set_text(r"Alt. [km]")
    ax1.axis["bottom"].label.set_text(r"Ground range [km]")
    ax1.invert_xaxis()

    ax1.minground   = minground
    ax1.maxground   = maxground
    ax1.minalt      = minalt
    ax1.maxalt      = maxalt
    ax1.Re          = Re

    fig.add_subplot(ax1, transform=tr)

    # create a parasite axes whose transData in RA, cz
    aux_ax          = ax1.get_aux_axes(tr)

    # for aux_ax to have a clip path as in ax
    aux_ax.patch    = ax1.patch

    # but this has a side effect that the patch is drawn twice, and possibly
    # over some other artists. So, we decrease the zorder a bit to prevent this.
    ax1.patch.zorder=0.9

    return ax1, aux_ax

def plot_rays(tx_lat,tx_lon,ranges,heights,
        maxground=None, maxalt=None,Re=6371,
        iono_arr=None,iono_param=None,
        iono_cmap='viridis', iono_lim=None, iono_title='Ionospheric Parameter',
        plot_rays=True,
        ray_path_data=None, 
        srch_ray_path_data=None, 
        fig=None, rect=111, ax=None, aax=None, cbax=None,
        plot_colorbar=True,
        iono_rasterize=False,scale_Re=1.,scale_heights=1.,terminator=False,
        title=None,**kwargs):
    """
    Plot a 2d ionospheric profile along a path.
    """
    
    Re_km       = Re
    heights_km  = heights
    Re          = scale_Re*Re
    heights     = scale_heights*heights

    if maxground is None:
        maxground = np.max(ranges)

    if maxalt is None:
        maxalt = np.max(heights)

    # Set up axes
    if not ax and not aax:
        ax, aax = curvedEarthAxes(fig=fig, rect=rect, 
            maxground=maxground, maxalt=maxalt,Re=Re,scale_heights=scale_heights)

    # Convert linear range into angular distance.
    thetas  = ranges/Re

    # Plot background ionosphere. ################################################## 
    if (iono_arr is not None) or (iono_param is not None):
        if iono_param == 'iono_en_grid' or iono_param == 'iono_en_grid_5':
            if iono_lim is None: iono_lim = (10,12)
            if iono_title == 'Ionospheric Parameter':
                iono_title = r"N$_{el}$ [$\log_{10}(m^{-3})$]"
            # Get the log10 and convert Ne from cm**(-3) to m**(-3)
            iono_arr    = np.log10(iono_arr*100**3)
            iono_arr[~np.isfinite(iono_arr)] = 0
        elif iono_param == 'iono_pf_grid' or iono_param == 'iono_pf_grid_5':
            if iono_lim is None: iono_lim = (0,10)
            if iono_title == 'Ionospheric Parameter':
                iono_title = "Plasma Frequency\n[MHz]"
        elif iono_param == 'collision_freq':
            if iono_lim is None: iono_lim = (0,8)
            if iono_title == 'Ionospheric Parameter':
                iono_title = r"$\nu$ [$\log_{10}(\mathrm{Hz})$]"
            iono_arr    = np.log10(iono_arr)

        if iono_lim is None:
            iono_mean   = np.mean(iono_arr)
            iono_std    = np.std(iono_arr)

            iono_0      = 50000
            iono_1      = 90000

            iono_lim    = (iono_0, iono_1)

        X, Y    = np.meshgrid(thetas,heights+Re)
        im      = aax.pcolormesh(X, Y, iono_arr[:-1,:-1],
                    vmin=iono_lim[0], vmax=iono_lim[1],zorder=1,
                    cmap=iono_cmap,rasterized=iono_rasterize,shading='auto')

        # Add a colorbar
        if plot_colorbar:
            cbar = plt.colorbar(im,label=iono_title,ax=ax,pad=0.02,
                    shrink=0.7,aspect=10)
            cbax = cbar.ax

    if terminator:
        # Terminator Resolution
#        term_dr = 50 # km
#        term_dh = 50 # km

        # Compute heights and ranges for terminator calculation
        term_heights    = heights_km.copy()
        term_shape      = (ranges.size,term_heights.size)

        # Convert ranges into lat/lons
        azm             = kwargs.get('azm')
        term_latlons    = geopack.greatCircleMove(tx_lat,tx_lon,ranges,azm)
        term_lats       = term_latlons[0]
        term_lons       = term_latlons[1]

        # Reshape heights and lat/lons into 2D arrays for array calculations.
        term_heights.shape  = (1,term_heights.size)
        term_lats.shape     = (term_lats.size,1)
        term_lons.shape     = (term_lons.size,1)

        TERM_HEIGHTS    = np.broadcast_to(term_heights,term_shape)
        TERM_LATS       = np.broadcast_to(term_lats,term_shape)
        TERM_LONS       = np.broadcast_to(term_lons,term_shape)

        # Use astropy to compute solar elevation angle at each point.
        obs     = ap_EarthLocation(lat=TERM_LATS*ap_u.deg,lon=TERM_LONS*ap_u.deg,
                    height=TERM_HEIGHTS*1e3*ap_u.m)

        ut      = ap_Time(kwargs.get('date'))
        frame   = ap_AltAz(obstime=ut,location=obs)
        az_el   = ap_get_sun(ut).transform_to(frame)
        el      = az_el.alt.value

        # Calculate horizon angle for each point assuming Re = 6371 km
        hzn     = np.degrees(np.arcsin(Re_km/(Re_km+TERM_HEIGHTS))) - 90.

        # Calculate delta horizon (number of degrees the sun is above the horizon)
        d_hzn   = el - hzn

        # Civil Twilight:         6 deg below horizon
        # Nautical Twilight:     12 deg below horizon
        # Astronomical Twilight: 18 deg below horizon

        twilight_deg    = 0.
        terminator_tf   = d_hzn <= twilight_deg
        terminator_tf   = (terminator_tf.astype(float)).T
        terminator_tf   = np.where(terminator_tf == 0, np.nan, terminator_tf)

        im      = aax.pcolormesh(X, Y, terminator_tf[:-1,:-1],
                    vmin=0, vmax=1., cmap='binary',zorder=2,alpha=0.10,
                    rasterized=iono_rasterize,shading='auto')

    # Plot Ray Paths ###############################################################
    if plot_rays:
        freq_s  = 'None'
        if ray_path_data is not None:
            rpd         = ray_path_data
            for inx,ray in enumerate(rpd):
                xx  = ray['ground_range']/Re
                yy  = ray['height']*scale_heights + Re
                aax.plot(xx,yy,color='white',lw=1.00,zorder=10)
            f   = ray['frequency']
            freq_s  = '{:0.3f} MHz'.format(float(f))

        if srch_ray_path_data is not None:
            rpd         = srch_ray_path_data
            for inx,ray in enumerate(rpd):
                xx  = ray['ground_range']/Re
                yy  = ray['height']*scale_heights + Re
                aax.plot(xx,yy,color='red',zorder=100,lw=7)

    # Plot Receiver ################################################################ 
    if 'rx_lat' in kwargs and 'rx_lon' in kwargs:
        rx_lat      = kwargs.get('rx_lat')
        rx_lon      = kwargs.get('rx_lon')
        rx_label    = kwargs.get('rx_label','Receiver')

        rx_theta    = kwargs.get('rx_range')/Re
        
        hndl    = aax.scatter([rx_theta],[Re],s=950,marker='*',color='red',ec='k',zorder=100,clip_on=False,label=rx_label)
        aax.legend([hndl],[rx_label],loc='upper right',scatterpoints=1,fontsize='large',labelcolor='black')
    
    if title is None:
        # Add titles and other important information.
        date_s      = kwargs.get('date').strftime('%Y %b %d %H:%M UT')
        tx_lat_s    = '{:0.2f}'.format(tx_lat) + r'$^{\circ}$N'
        tx_lon_s    = '{:0.2f}'.format(tx_lon) + r'$^{\circ}$E'
        azm_s       = '{:0.1f}'.format(kwargs['azm'])   + r'$^{\circ}$'
        if plot_rays:
            ax.set_title('{} - {}'.format(date_s,freq_s),fontsize=title_fontsize)
        else:
            ax.set_title('{}'.format(date_s),fontsize=title_fontsize)
    else:
        ax.set_title(title,fontsize=title_fontsize)

    return ax, aax, cbax

class RayTraceAndPlot(object):
    def __init__(self,iono_nc=None, freq=14, nhops=1):
        if iono_nc is None:
        #    iono_nc = 'data/iri_tid_300km/20181012.1830-20181012.1830_FHE__profile.nc'
            iono_nc = 'data/iri_tid_1000km/20181512.2000-20181512.2000_TX__profile.nc'

        print('Loading ionospheric grid {!s}... '.format(iono_nc))
        iono_ds         = xr.load_dataset(iono_nc)
        self.iono_ds    = iono_ds

        self.prmd = prmd = {}
        prmd['UT']             = pd.to_datetime(iono_ds['date'].values[0])
        prmd['origin_lat']     = iono_ds.attrs['tx_lat']
        prmd['origin_lon']     = iono_ds.attrs['tx_lon']
        prmd['ray_bear']       = np.round(iono_ds.attrs['azm'],1)
        prmd['start_height']   = float(iono_ds['alt'].min())
        prmd['height_inc']     = np.diff(iono_ds['alt'])[0]
        prmd['range_inc']      = np.diff(iono_ds['range'])[0]
        prmd['heights']        = iono_ds['alt'].values
        prmd['ranges']         = iono_ds['range'].values

        en_ds                  = np.squeeze(iono_ds['electron_density'].values).T # Needs to be dimensions of [height, range]
        en_ds[en_ds < 0 ]      = 0                                                # Make sure all electron densities are >= 0.
        en_ds                  = en_ds / (100**3)                                 # PyLap needs electron densities in electrons per cubic cm.
        prmd['iono_en_grid']   = en_ds
        prmd['iono_en_grid_5'] = en_ds                                            # We are not calculating Doppler shift, so the 5 minute electron densities can be the same as iono_en_grid
        prmd['collision_freq'] = en_ds*0.                                         # Ignoring collision frequencies
        prmd['irreg']          = np.zeros([4,en_ds.shape[1]])                     # Ignoring ionospheric irregularities.

        prmd['elevs']          = np.arange(2, 62, 0.5, dtype = float)             # py
        prmd['num_elevs']      = len(prmd['elevs'])
        prmd['freq']           = freq                                             # Ray Frequency (MHz)
        prmd['freqs']          = freq * np.ones(prmd['num_elevs'], dtype = float) # Need to pass a vector of frequencies the same length as len(elevs)
        prmd['tol']            = [1e-7, 0.01, 10]                                 # ODE tolerance and min/max step sizes
        prmd['nhops']          = nhops                                            # number of hops to raytrace
        prmd['irregs_flag']    = 0                                                # no irregularities - not interested in Doppler spread or field aligned irregularities

        self.ray_trace()

    def ray_trace(self):
        print('Generating {} 2D NRT rays ...'.format(self.prmd['num_elevs']))
        rt_keys = ['origin_lat', 'origin_lon', 'elevs', 'ray_bear', 'freqs', 'nhops', 'tol', 'irregs_flag', 'iono_en_grid', 'iono_en_grid_5', 'collision_freq', 'start_height', 'height_inc', 'range_inc', 'irreg']
        rt_prms = [self.prmd[key] for key in rt_keys]
        ray_data, ray_path_data, ray_path_state = raytrace_2d(*rt_prms)

        self.ray_data       = ray_data
        self.ray_path_data  = ray_path_data
        self.ray_path_state = ray_path_state

    def plot_figure(self,fpath='output.png',figsize=(40,10),**kwargs):
        fig = plt.figure(figsize=figsize)
        self.plot_ax(fig=fig,**kwargs)

        print('Saving Figure: {!s}'.format(fpath))
        fig.savefig(fpath,bbox_inches='tight')

        # Remove whitespace using mogrify since the curved axes are
        # not compatible with bbox_inches='tight'
        try:
            cmd = f'mogrify -trim {fpath}'
            os.system(cmd)
        except:
            print(f'ERROR: Could not run {cmd}')


    def plot_ax(self,fig=None,end_range=3000,end_ht=500,
            cbax                = None,
            plot_colorbar       = True,
            panel_rect          = None,
            title_size          = None,
            ticklabel_size      = None,
            label_size          = None,
            cbar_ticklabel_size = None,
            cbar_label_size     = None,
            **kwargs):

        if ticklabel_size is not None:
            mpl.rcParams['xtick.labelsize']         = ticklabel_size
            mpl.rcParams['ytick.labelsize']         = ticklabel_size

        prmd    = self.prmd

        _pltd    = {}
        _pltd['tx_lat']             = prmd['origin_lat']
        _pltd['tx_lon']             = prmd['origin_lon']
        _pltd['ranges']             = prmd['ranges']
        _pltd['heights']            = prmd['heights']
        _pltd['maxground']          = end_range
        _pltd['maxalt']             = end_ht
        _pltd['date']               = prmd['UT']    # Used for terminator.
        _pltd['azm']                = prmd['ray_bear']# Used for terminator.
        _pltd['iono_arr']           = prmd['iono_en_grid']
        _pltd['iono_param']         = 'iono_en_grid'
        _pltd['ray_path_data']      = self.ray_path_data
        _pltd['fig']                = fig 
        _pltd['cbax']               = cbax
        _pltd['plot_colorbar']      = plot_colorbar
        _pltd['title']              = ''
        _pltd.update(kwargs)

        ax, aax, cbax   = plot_rays(**_pltd)

        title_prms = {'weight':'bold'}
        if title_size is not None:
            title_prms.update({'size':title_size})

        title   = []
        tid_str = ' ($\lambda_h$=1000 km)'
        title.append('IRI2016 Perturbed with LSTID' + tid_str)
        title.append('{!s}'.format(prmd['UT'].strftime('%Y %b %d %H:%M UTC')))
        title   = '\n'.join(title)
        ax.set_title(title,loc='left',**title_prms)

        title   = []
        title.append('{!s} MHz Raytrace'.format(prmd['freq']))
        title.append('Origin {:0.1f}\N{DEGREE SIGN}N, {:0.1f}\N{DEGREE SIGN}E, {:0.0f}\N{DEGREE SIGN} AZM'.format(prmd['origin_lat'],prmd['origin_lon'],prmd['ray_bear']))
        title   = '\n'.join(title)
        ax.set_title(title,loc='right',**title_prms)

        if label_size is not None:
            ax.xaxis.label.set_size(label_size)
            ax.yaxis.label.set_size(label_size)

        if cbar_ticklabel_size is not None:
            for ytl in cbax.get_yticklabels():
                ytl.set_size(cbar_ticklabel_size)

        if cbar_label_size is not None:
            cbax.xaxis.label.set_size(cbar_label_size)
            cbax.yaxis.label.set_size(cbar_label_size)

        if panel_rect is not None:
            p_x00   = panel_rect[0] + 0.025
            p_y00   = panel_rect[1] - 0.01
            p_wd    = panel_rect[2]
            p_ht    = panel_rect[3]

            rt_wd   = 0.85      * p_wd
            rt_x00  = p_x00

            cb_wd   = (1-rt_wd) * p_wd
            cb_x00  = rt_x00 + rt_wd + 0.025

            rt_rect = [rt_x00, p_y00, rt_wd, p_ht]
            cb_rect = [cb_x00, p_y00, cb_wd, p_ht]
            ax.set_position(rt_rect)
            cbax.set_position(cb_rect)

        result = {}
        result['ax']    = ax
        result['aax']   = aax
        result['cbax']  = cbax
        return result

if __name__ == '__main__':
    iono_nc = 'data/iri_tid_1000km/20181512.2000-20181512.2000_TX__profile.nc'
    output_dir  = os.path.join('output','Figure_2b')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fname   = 'Figure_2b_TID_RayTrace.png'
    fpath   = os.path.join(output_dir,fname)

    RTaP    = RayTraceAndPlot(iono_nc)
    RTaP.plot_figure(fpath=fpath)
