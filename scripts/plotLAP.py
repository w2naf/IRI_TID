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

# The GeographicLib is likely more accurate than geooack.
# Geographiclib - https://geographiclib.sourceforge.io/Python/2.0/
# mamba install conda-forge::geographiclib
from geographiclib.geodesic import Geodesic
geod = Geodesic.WGS84

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
            cbar = plt.colorbar(im,label=iono_title,ax=ax,pad=0.02)
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
                aax.plot(xx,yy,color='red',zorder=100,lw=4)

    # Plot Receiver ################################################################ 
    if 'rx_lat' in kwargs and 'rx_lon' in kwargs:
        rx_lat      = kwargs.get('rx_lat')
        rx_lon      = kwargs.get('rx_lon')
        rx_label    = kwargs.get('rx_label','Receiver')

        # Determine the ranges and azimuth along the profile path.
        invl    = geod.InverseLine(tx_lat,tx_lon,rx_lat,rx_lon)
        rx_dist_km = invl.s13*1e-3   # Distance in km
        rx_az      = invl.azi1

        rx_theta    = rx_dist_km/Re
        
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

def main(iono_ds,freq=14.,
        elev_0 = 2, elev_1 = 62, delev = 0.5,
        plot_end_range  = 2000.,
        plot_end_height = 500.,
        png_fpath='raytrace.png'):
    """
    freq:   Ray Frequency (MHz)
    elev_0: Minimum Elevation Angle (degrees)
    elev_1: Maximum Elevation Angle (degrees)
    delev:  Elevation Angle Step Size (degrees)
    plot_end_range:  End range of plot [km]
    plot_end_height: End height of plot [km]
    """
    UT              = pd.to_datetime(iono_ds['date'].values[0])
    origin_lat      = iono_ds.attrs['tx_lat']
    origin_lon      = iono_ds.attrs['tx_lon']
    ray_bear        = np.round(iono_ds.attrs['azm'],1)
    start_height    = float(iono_ds['alt'].min())
    height_inc      = np.diff(iono_ds['alt'])[0]
    range_inc       = np.diff(iono_ds['range'])[0]
    heights         = iono_ds['alt'].values
    ranges          = iono_ds['range'].values

    UT              = pd.to_datetime(iono_ds['date'].values[0])
    origin_lat      = iono_ds.attrs['tx_lat']
    origin_lon      = iono_ds.attrs['tx_lon']
    ray_bear        = np.round(iono_ds.attrs['azm'],1)
    start_height    = float(iono_ds['alt'].min())
    height_inc      = np.diff(iono_ds['alt'])[0]
    range_inc       = np.diff(iono_ds['range'])[0]
    heights         = iono_ds['alt'].values
    ranges          = iono_ds['range'].values

    en_ds           = np.squeeze(iono_ds['electron_density'].values).T  #Needs to be dimensions of [height, range]
    en_ds[en_ds < 0 ] = 0 # Make sure all electron densities are >= 0.
    en_ds           = en_ds / (100**3) # PyLap needs electron densities in electrons per cubic cm.
    iono_en_grid    = en_ds
    iono_en_grid_5  = iono_en_grid      # We are not calculating Doppler shift, so the 5 minute electron densities can be the same as iono_en_grid
    collision_freq  = iono_en_grid*0.   # Ignoring collision frequencies
    irreg           = np.zeros([4,en_ds.shape[1]])  # Ignoring ionospheric irregularities.

    elevs           = np.arange(elev_0, elev_1, delev, dtype = float) # py
    num_elevs       = len(elevs)
    freqs           = freq * np.ones(num_elevs, dtype = float) # Need to pass a vector of frequencies the same length as len(elevs)
    tol             = [1e-7, 0.01, 10]  # ODE tolerance and min/max step sizes
    nhops           = 2                 # number of hops to raytrace
    irregs_flag     = 0                 # no irregularities - not interested in Doppler spread or field aligned irregularities

#    import ipdb; ipdb.set_trace()
    print('Generating {} 2D NRT rays ...'.format(num_elevs))
    ray_data, ray_path_data, ray_path_state = \
       raytrace_2d(origin_lat, origin_lon, elevs, ray_bear, freqs, nhops,
                   tol, irregs_flag, iono_en_grid, iono_en_grid_5,
               collision_freq, start_height, height_inc, range_inc, irreg)

    tx_lat      = iono_ds.attrs.get('tx_lat')
    tx_lon      = iono_ds.attrs.get('tx_lon')
    tx_call     = iono_ds.attrs.get('tx_call')
    rx_lat      = iono_ds.attrs.get('rx_lat')
    rx_lon      = iono_ds.attrs.get('rx_lon')
    rx_call     = iono_ds.attrs.get('rx_call')

    ########################
    # Convert to DataFrame #
    ########################
    rd_df = {}
    rd_df['ray_id'] = []
    rd_df['hop_id'] = []
    keys    = list(ray_data[0].keys())
    for key in keys:
        rd_df[key] = []

    for ray_id,rd in enumerate(ray_data):
        nhops   = len(rd['ray_label'])

        rd_df['ray_id'] += [ray_id]*nhops
        rd_df['hop_id'] += list(range(nhops))

        for key in keys:
            if key in ['frequency','nhops_attempted']:
                val = [rd[key]] * nhops
            else:
                val = rd[key].tolist()

            rd_df[key] += val

    rd_df   	= pd.DataFrame(rd_df)
    csv_fpath	= png_fpath[:-4]+'.rayData.csv'

    hdr = []
    hdr.append('# '+csv_fpath)
    hdr.append('# TX: {!s} {:0.1f}\N{DEGREE SIGN}N, {:0.1f}\N{DEGREE SIGN}E, {:0.0f}\N{DEGREE SIGN} AZM'.format(tx_call,origin_lat,origin_lon,ray_bear))
	
    line= """#
# Ray Data CSV File Produced from PyLAP/PHaRLAP 2D Ray Tracing
#
# lat: Geographic latitude of Ray Trace Termination
# lon: Geographic longitude of Ray Trace Termination
# ground_range: [km] Ground range is the distance along the ground from the ray transmit location to the point on the ground where the ray lands.
# group_range: [km] Group range is the time taken for the ray to traverse the path multiplied by c, the speed of light.
#       Due to the group retardation effect of the ionosphere the actual speed of the wave packet is < c and so the
#       group range will be greater than the actual physical path of the ray.
# phase_path: [km] Phase path is the accumulated phase over the path converted to a distance in free space. The phase speed is > c and so
#       the phase path will be less than actual physical path of the ray.
# geometric_path_length: [km] Geometric path length is the actual distance travelled by the radio waves
# initial_elev: [degrees] Take-off elevation angle
# final_elev: [degrees] Received elevation angle
# apogee: [km] Altitude of highest point of the ray.
# gnd_rng_to_apogee: [km] Ground distance from transmitter to apogee location.
# plasma_frequency_at_apogee: [MHz]
# virtual_height: [km] Virtual Height at Apogee (assumes specular reflection)
# effective_range [m] Effective range accounts for the focussing/defocussing effects.
#       Davies, 1990, “Ionospheric Radio” gives a good overview of group path and phase path in chapter 1. Effective range is
#       discussed in Chapter 7.3.2 (page 210 – 213)
# total_absorption:
# deviative_absorption:
# TEC_path: 
# Doppler_shift:
# Doppler_spread: 
# FAI_backscatter_loss:
# frequency: [MHz] Frequency of radio wave being raytraced.
# nhops_attempted:
# ray_label: - label for each hop attempted which indicates what the ray has done.
#        1  for ray reaching ground
#        0  for ray becoming evanescent, raytracing terminated
#       -1  for field aligned backscatter - ray reflected with appropriate scattering loss, raytracing terminated
#       -2  ray has penetrated the ionosphere - raytracing terminated
#       -3  ray has exceeded max. ground range - raytracing terminated
#       -4  ray angular coordinate has become negative (bad - should never happen) - raytracing terminated
#       -5  ray has exceeded the maximum allowed points along path (20000 points) - raytracing terminated
#       -6  ray is near antipodal point, the WGS84 coordinate conversion routines are unreliable - terminate raytracing
#     -100  a catastrophic error occured - terminate raytracing
#
"""
    hdr.append(line)
    header = '\n'.join(hdr)
    with open(csv_fpath,'w') as fl:
        fl.write(header)
    rd_df.to_csv(csv_fpath,mode='a',index=False)

    ########################
    # Search for Rays      #
    ########################

    # Determine the ranges and azimuth along the profile path.
    invl        = geod.InverseLine(tx_lat,tx_lon,rx_lat,rx_lon)
    rx_dist_km = invl.s13*1e-3   # Distance in km

    search_radius       = 20.
    srch_ray_path_data  = None
    if search_radius is not None:
        srch_range_0    = rx_dist_km - search_radius
        srch_range_1    = rx_dist_km + search_radius

        tf_0    = rd_df['ground_range'] >= srch_range_0
        tf_1    = rd_df['ground_range'] <  srch_range_1
        tf_2    = rd_df['ray_label'] == 1
        tf      = np.logical_and.reduce([tf_0,tf_1,tf_2])

        srch_rd_df  = rd_df[tf]

        if len(srch_rd_df) > 0:
            srch_ray_inx    = srch_rd_df['ray_id'].unique()
            srch_ray_path_data = [ray_path_data[x] for x in srch_ray_inx]

            hdr.append('# RX: {!s} {:0.1f}\N{DEGREE SIGN}N, {:0.1f}\N{DEGREE SIGN}E'.format(rx_call,rx_lat,rx_lon))
            hdr.append('# Search Radius: {:0.0f} km'.format(search_radius))
            header = '\n'.join(hdr)

            srch_csv_fpath	= png_fpath[:-4]+'.search_rayData.csv'
            with open(srch_csv_fpath,'w') as fl:
                fl.write(header)
            srch_rd_df.to_csv(srch_csv_fpath,mode='a',index=False)

    ###################
    ### Plot Result ###
    ###################

    kwRays = {}
    kwRays['rx_lat']    = iono_ds.attrs.get('rx_lat')
    kwRays['rx_lon']    = iono_ds.attrs.get('rx_lon')
    kwRays['rx_label']  = iono_ds.attrs.get('rx_call')

    fig = plt.figure(figsize=(35,10))
    ax, aax, cbax   = plot_rays(origin_lat,origin_lon,ranges,heights,
            maxground=plot_end_range, maxalt=plot_end_height,Re=6371,date=UT,azm=ray_bear,
            iono_arr=iono_en_grid,iono_param='iono_en_grid',
            iono_cmap='viridis', iono_lim=None, iono_title='Ionospheric Parameter',
            plot_rays=True,
            ray_path_data=ray_path_data, 
            srch_ray_path_data=srch_ray_path_data, 
            fig=None, rect=111, ax=None, aax=None, cbax=None,
            plot_colorbar=True,title='',
            iono_rasterize=False,scale_Re=1.,scale_heights=1.,terminator=False,**kwRays)

    title   = []
    title.append(os.path.basename(png_fpath))
    title.append('{!s}'.format(UT.strftime('%Y %b %d %H:%M UTC')))
    title   = '\n'.join(title)
    ax.set_title(title,loc='left')

    title   = []
    title.append('TX: {!s} {:0.1f}\N{DEGREE SIGN}N, {:0.1f}\N{DEGREE SIGN}E, {:0.0f}\N{DEGREE SIGN} AZM'.format(tx_call,origin_lat,origin_lon,ray_bear))
    title.append('RX: {!s} {:0.1f}\N{DEGREE SIGN}N, {:0.1f}\N{DEGREE SIGN}E'.format(rx_call,rx_lat,rx_lon))
    if search_radius is not None:
        title.append('Search Radius: {:0.0f} km'.format(search_radius))
    title   = '\n'.join(title)
    ax.set_title(title,loc='right')

    print('Saving Figure: {!s}'.format(png_fpath))
    fig.savefig(png_fpath,bbox_inches='tight')

    return ray_data, ray_path_data, ray_path_state
