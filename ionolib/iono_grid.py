#!/usr/bin/env python
import os
import shutil
import datetime
import bz2

from multiprocessing import Pool

import tqdm

import matplotlib

from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon

from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import scipy.interpolate
import pandas as pd
import xarray as xr

import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# The GeographicLib is likely more accurate than geooack.
# Geographiclib - https://geographiclib.sourceforge.io/Python/2.0/
# mamba install conda-forge::geographiclib
from geographiclib.geodesic import Geodesic
geod = Geodesic.WGS84

# Geopack is vectorized and faster than geographiclib, so it is still useful
# for certain applications.
try:
    from . import geopack
except:
    import geopack

Re = 6371 # Earth Radius in km

try:
    # Michael Hirsch's IRI2016 - https://github.com/space-physics/iri2016 
    import iri2016
except:
    print('Michael Hirsch\'s IRI2016 Not Installed')
    print('If you want to use this engine, please install from:')
    print('https://github.com/space-physics/iri2016') 

try:
    # Victoria Forsythe's PyIRI - https://github.com/victoriyaforsythe/PyIRI
    import PyIRI
    import PyIRI.main_library as ml
except:
    print('Victoria Forsythe\'s PyIRI Not Installed')
    print('If you want to use this engine, please install from:')
    print('https://github.com/victoriyaforsythe/PyIRI')


def calculate_scale(data,stddevs=2.,lim='auto'):
    if lim == 'auto':
        mean    = np.nanmean(np.abs(data))
        std     = np.nanstd(np.abs(data))
        lim     = mean + stddevs*std
        scale   = (0,lim)
    else:
        scale   = lim

    ticks   = None
    cmap    = matplotlib.cm.viridis

    return scale,ticks,cmap

class iono_3d(object):
    def __init__(self,
            sDate = datetime.datetime(2016,6,21),
            eDate = datetime.datetime(2016,6,21,3),
            dt    = datetime.timedelta(hours=1),
            hgt_0 =  250., hgt_1 =  300., hgt_step=100.0,
            lat_0 =   30., lat_1 =   80., lat_step=5.0,
            lon_0 = -100., lon_1 =  -40., lon_step=5.0,
            engine = "PyIRI",
            data_dir = 'data', cache = True):
        """
        engine:
            'PyIRI':    Victoria Forsythe's PyIRI (https://github.com/victoriyaforsythe/PyIRI)
            'iri2016':  Michael Hirsch's IRI2016 Python Wrapper (https://github.com/space-physics/iri2016)
        """

        dates   = [sDate]
        while dates[-1] < eDate:
            dates.append(dates[-1]+dt)

        lats    = np.arange(lat_0,lat_1,lat_step)
        lons    = np.arange(lon_0,lon_1,lon_step)
        alts    = np.arange(hgt_0,hgt_1,hgt_step)

        self.engine     = engine
        self.dates      = dates
        self.lats       = lats
#        self.lons       = lons % 360. # Adjust so lons are between 0 and 360 deg.
        self.lons       = lons
        self.alts       = alts
        self.lat_step   = lat_step
        self.lon_step   = lon_step
        self.alt_step   = hgt_step

        fname   = []
        fname.append(engine)
        fname.append(sDate.strftime('%Y%m%d.%H%MUT'))
        fname.append(eDate.strftime('%Y%m%d.%H%MUT'))
        fname.append('{:.1f}min'.format(dt.total_seconds()/60.))
        fname.append('{:.0f}km'.format(hgt_0))
        fname.append('{:.0f}km'.format(hgt_1))
        fname.append('{:.0f}km'.format(hgt_step))
        fname.append('{:.0f}N'.format(lat_0))
        fname.append('{:.0f}N'.format(lat_1))
        fname.append('{:.0f}deg'.format(lat_step))
        fname.append('{:.0f}E'.format(lon_0))
        fname.append('{:.0f}E'.format(lon_1))
        fname.append('{:.0f}deg'.format(lon_step))
        fname       = '_'.join(fname)+'.nc'
        self.fname  = fname
        fpath       = os.path.join(data_dir,fname)

        if cache:
            if not os.path.exists(data_dir):
                os.makedirs(data_dir)

            if os.path.exists(fpath):
                self.iri_dataset = xr.load_dataset(fpath)
            else:
                ds  = self.run_iri() 
                ds.to_netcdf(fpath)
        else:
            ds  = self.run_iri() 

        self.profiles   = {}

    def run_iri(self):
        """
        Generate 3D array of electron densities from IRI.
        """

        nDates  = len(self.dates)
        nLats   = len(self.lats)
        nLons   = len(self.lons)
        nAlts   = len(self.alts)

        edens   = np.zeros([nDates,nLats,nLons,nAlts])

        ################################################################################
        iri_rundcts = []
        edens_inxs  = []
        for dInx,date in tqdm.tqdm(enumerate(self.dates),desc='Prepping IRI Run Dictionaries',dynamic_ncols=True,total=len(self.dates)):
            if self.engine == 'PyIRI':
                year    = date.year
                month   = date.month
                day     = date.day
                dhour   = date.hour + date.minute/60.
                ahr     = np.array([dhour])
                alat_2d, alon_2d = np.meshgrid(self.lats,self.lons)
                alat    = np.reshape(alat_2d, alat_2d.size)
                alon    = np.reshape(alon_2d, alon_2d.size)
                aalt    = self.alts
                f107    = 100

                # Specify what coefficients to use for the peak of F2 layer:
                # 0 = CCIR, 1 = URSI
                ccir_or_ursi    = 0

                f2, f1, e_peak, es_peak, sun, mag, edp = ml.IRI_density_1day(year, month, day, 
                        ahr, alon, alat, aalt, f107, PyIRI.coeff_dir, ccir_or_ursi)

                for inx,(lat,lon) in enumerate(zip(alat,alon)):
                    latinx  = np.where(self.lats == lat)[0][0]
                    loninx  = np.where(self.lons == lon)[0][0]
                    edens[dInx,latinx,loninx,:]  = edp[0,:,inx]

            elif self.engine == 'iri2016':
                for latinx, lat in enumerate(self.lats):
                    for loninx, lon in enumerate(self.lons):
                        
                        altkmrange  = [np.min(self.alts),np.max(self.alts),self.alt_step]
                        iri = iri2016.IRI(date, altkmrange, lat, lon)
                        # edens[dinx,latinx,loninx,alt_inx]
                        edens[dInx,latinx,loninx,:] = iri['ne']

        ################################################################################

        # Create XArray Dataset
        ds  = xr.Dataset(
                data_vars=dict(
                    electron_density=(['date','lat','lon','alt'],edens)
                    ),
                coords=dict(
                    date    = self.dates,
                    lat     = self.lats,
                    lon     = self.lons,
                    alt     = self.alts
                    ),
                attrs=dict(description='IRI Run')
                )
        self.iri_dataset    = ds
        return ds

    def generate_tx_rx_profile(self,tx_lat,tx_lon,rx_lat,rx_lon,range_step=5,max_range=None,
            tx_call='None',rx_call='None',interp_type='nearest',attrs={}):
        """
        Creates a 2D slice profile of the 3D IRI grid between a transmit location and receive location for each
        time step in the grid.
        tx_lat:         Transmitter (starting) latitude
        tx_lon:         Transmitter (starting) longitude
        rx_lat:         Receiver (ending) latitude
        rx_lon:         Receiver (ending) longitude
        range_step:     Range step in km
        max_range:      Override distance between TX and RX if not None. [km]
        tx_call:        String identifying the transmitter
        rx_call:        String identifying the receiver
        interp_type:    Method of interpolation to use.
                            'nearest':  Nearest neighbor (fastest)
                            'linear':   Linear (smoothest)
        """
        # Create meshgrid of original coordinates.
        LATS, LONS, ALTS    = np.meshgrid(self.lats,self.lons,self.alts,indexing='ij')
        coords              = list(zip(LATS.flatten(),LONS.flatten(),ALTS.flatten()))

        # Determine the ranges and azimuth along the profile path.
        invl    = geod.InverseLine(tx_lat,tx_lon,rx_lat,rx_lon)
        dist    = invl.s13*1e-3   # Distance in km
        az      = invl.azi1

        if max_range is not None:
            dist = max_range

        ranges  = np.arange(0,dist,range_step)

        glats   = []
        glons   = []
        for x in ranges:
            s   = min(x*1e3,invl.s13) # invl.s13 is the total line distance in m
            tmp = invl.Position(s,Geodesic.STANDARD)
            glat        = tmp['lat2']
            glon        = tmp['lon2']

            glats.append(glat)
            glons.append(glon)
        flats_flons = np.array([glats,glons]).T

        # Put the field points into a mesh.
        fLATS   = np.zeros([len(ranges),len(self.alts)])
        fLONS   = np.zeros([len(ranges),len(self.alts)])
        fALTS   = np.zeros([len(ranges),len(self.alts)])
        for rInx,flat_flon in enumerate(flats_flons):
            fLATS[rInx,:]   = flat_flon[0]
            fLONS[rInx,:]   = flat_flon[1]
            fALTS[rInx,:]   = self.alts

        dict_key    = '{}-{}'.format(tx_call,rx_call)
        Ne_profile  = np.zeros([len(self.dates),len(ranges),len(self.alts)])
        for dateInx,date in tqdm.tqdm(enumerate(self.dates),desc='Interpolating Profiles',dynamic_ncols=True):
            tqdm.tqdm.write('INTERP: {!s}'.format(date))
            this_edens          = self.iri_dataset['electron_density'].values[dateInx,:,:,:]

            if interp_type == 'nearest':
                interp  = sp.interpolate.NearestNDInterpolator(coords,this_edens.flatten())
            elif interp_type == 'linear':
                interp  = sp.interpolate.LinearNDInterpolator(coords,this_edens.flatten())

            edens_profile       = interp(fLATS,fLONS,fALTS)
            Ne_profile[dateInx,:,:]    = edens_profile

        # Calculate range from start as an angle [radians]
        # Computed the same way as in raydarn fortran code.
        field_lats  = flats_flons[:,0]
        field_lons  = flats_flons[:,1]
        edensTHT    = np.arccos( np.cos(field_lats[0]*np.pi/180.)*np.cos(field_lats*np.pi/180.)* \
                            np.cos((field_lons - field_lons[0])*np.pi/180.) \
                    + np.sin(field_lats[0]*np.pi/180.)*np.sin(field_lats*np.pi/180.))

        # Set initial edensTHT = 0
        edensTHT[0] = 0

        # Base filename to be used with this profile.
        fname_base  = '{tx_call}_{rx_call}'.format(tx_call=tx_call,rx_call=rx_call)

        _attrs = dict(
                    tx_call     = tx_call,
                    tx_lat      = tx_lat,
                    tx_lon      = tx_lon,
                    rx_call     = rx_call,
                    rx_lat      = rx_lat,
                    rx_lon      = rx_lon,
                    azm         = az,
                    fname_base  = fname_base
                )

        _attrs.update(attrs)

        # Create XArray Dataset
        ds  = xr.Dataset(
                data_vars=dict(
                    electron_density = (['date','range','alt'],Ne_profile),
                    dip              = (['date','range','alt'],Ne_profile*0),
                    glats            = (['range'],field_lats),
                    glons            = (['range'],field_lons),
                    edensTHT         = (['range'],edensTHT)
                    ),
                coords=dict(
                    date        = self.dates,
                    range       = ranges,
                    alt         = self.alts
                    ),
                attrs=_attrs
            )

        self.profiles[dict_key] = ds
        return self.profiles

    def profiles_to_netcdf(self,keys=None,output_dir=None):
        if keys is None:
            keys = self.profiles.keys()
        
        paths = []
        for key in keys:
            profl       = self.profiles[key]
            dates       = list(map(pd.to_datetime,profl['date'].values))
            fname_base  = profl.attrs['fname_base']

            dminS = min(dates).strftime('%Y%d%m.%H%M')
            dmaxS = max(dates).strftime('%Y%d%m.%H%M')
            _filename = os.path.join(output_dir,'{!s}-{!s}_{!s}_profile.nc'.format(dminS,dmaxS,fname_base))

            profl.to_netcdf(_filename)

    def plot_profiles(self,keys=None,output_dir='output',filename=None,figsize=(10,6)):

        if keys is None:
            keys = self.profiles.keys()
        
        paths = []
        for key in keys:
            profl       = self.profiles[key]
            edens       = profl['electron_density'].values
            dates       = list(map(pd.to_datetime,profl['date'].values))
            alts        = profl['alt'].values
            ranges      = profl['range'].values
            lats        = profl['glats']
            lons        = profl['glons']
            fname_base  = profl.attrs['fname_base']

            range_step  = np.mean(np.diff(ranges))
            alt_step    = np.mean(np.diff(alts))

            scale,ticks,cmap = calculate_scale(edens)
            bounds  = np.linspace(scale[0],scale[1],256)
            norm    = matplotlib.colors.BoundaryNorm(bounds,cmap.N)

            for dInx,date in tqdm.tqdm(enumerate(dates),desc='Plotting Profiles',dynamic_ncols=True):
                tqdm.tqdm.write('Plotting {!s}: {!s}'.format(key,date))
                this_edens  = edens[dInx,:,:]
                
                fig  = plt.figure(figsize=figsize)
                ax   = fig.add_subplot(111)

                pcoll = ax.pcolorfast(ranges,alts,this_edens.T,cmap=cmap,norm=norm)

                ax.set_xlim(ranges.min(), ranges.max())
                ax.set_ylim(alts.min(), alts.max())

                ax.set_xlabel('Range [km]')
                ax.set_ylabel('Altitude [km]')

                cbar    = fig.colorbar(pcoll,orientation='vertical',shrink=0.60,pad=.10,ticks=ticks)
                txt     = r'IRI Electron Density [m$^{-3}$]'
                cbar.set_label(txt)

                txt = []
                txt.append('IRI Electron Density')
                txt.append('{0} {1}'.format(key,date.strftime('%d %b %Y %H%M UT')))
                ax.set_title('\n'.join(txt))

                fig.tight_layout()

                if filename is None:
                    _filename = os.path.join(output_dir,'{!s}_{!s}_profile.png'.format(date,fname_base))
                else:
                    _filename = os.path.join(output_dir,filename)

                fig.savefig(_filename,bbox_inches='tight')
                plt.close()
                paths.append(_filename)
        return paths

    def plot_maps(self,alt=250.,output_dir='output',figsize=(15,8),
            xlim=None,ylim=None,plot_profile_paths='all'):
        ds      = self.iri_dataset
        edens   = ds['electron_density'].values
        lats    = ds['lat'].values
        lons    = ds['lon'].values
        alts    = ds['alt'].values
        dates   = list(map(pd.to_datetime,ds['date'].values))

        alt_inx = np.argmin(np.abs(alts-alt))

        scale,ticks,cmap = calculate_scale(edens[:,:,:,alt_inx])
        bounds      = np.linspace(scale[0],scale[1],256)
        norm        = matplotlib.colors.BoundaryNorm(bounds,cmap.N)

        for dInx,date in tqdm.tqdm(enumerate(dates),desc='Plotting Maps',dynamic_ncols=True):
            tqdm.tqdm.write('Plotting: {!s}'.format(date))
            this_edens  = edens[dInx,:,:]

            fig     = plt.figure(figsize=figsize)
            ax      = fig.add_subplot(111,projection=ccrs.PlateCarree())

    #        ax.coastlines(zorder=10,color='k')
    #        ax.add_feature(cartopy.feature.LAND)
    #        ax.add_feature(cartopy.feature.OCEAN)
    #        ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    #        ax.add_feature(cartopy.feature.RIVERS)
            ax.add_feature(cartopy.feature.COASTLINE)
            ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
            ax.set_title('')
            ax.gridlines(draw_labels=True)

            LONS, LATS  = np.meshgrid(lons,lats)
            data        = this_edens[:,:,alt_inx]
            pcoll       = ax.pcolormesh(LONS,LATS,data,cmap=cmap,norm=norm)

            if plot_profile_paths == 'all':
                for prof_key,profile in self.profiles.items():
                    glons       = profile['glons']
                    tf          = glons > 180
                    glons[tf]   = glons[tf]-360.
                    glats       = profile['glats']
                    ax.plot(glons,glats,marker='o',color='k',lw=4,zorder=100)

                    tx_call     = profile.attrs['tx_call']
                    tx_lat      = profile.attrs['tx_lat']
                    tx_lon      = profile.attrs['tx_lon']
                    if tx_lon > 180:
                        tx_lon = tx_lon - 360.

                    rx_call     = profile.attrs['rx_call']
                    rx_lat      = profile.attrs['rx_lat']
                    rx_lon      = profile.attrs['rx_lon']
                    if rx_lon > 180:
                        rx_lon = rx_lon - 360.

                    ax.scatter(tx_lon,tx_lat,marker='*',s=450,zorder=110,label=tx_call,ec='k',fc='red')
                    ax.scatter(rx_lon,rx_lat,marker='*',s=450,zorder=110,label=rx_call,ec='k',fc='orange')
                    fontdict = {'size':'x-large','weight':'bold'}
                    offset = 1.1
                    ax.text(tx_lon,tx_lat+offset,tx_call,fontdict=fontdict,ha='center')
                    ax.text(rx_lon,rx_lat+offset,rx_call,fontdict=fontdict,ha='center')

            # Plot Wave Source Locations
            if hasattr(self,'wave_list'):
                for wave in self.wave_list:
                    wlat    = wave['src_lat']
                    wlon    = wave['src_lon']
                    ax.scatter(wlon,wlat,marker='^',s=100,zorder=110,ec='k',fc='red')

            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            cbar_label  = r'IRI Electron Density [m$^{-3}$]'
            cbar        = fig.colorbar(pcoll,orientation='vertical',shrink=0.65,pad=0.075,ticks=ticks)
            cbar.set_label(cbar_label,fontdict={'weight':'bold','size':'large'})

            # Plot Title
            txt = []
            txt.append('{0} - Alt: {1:.0f} km'.format(date.strftime('%d %b %Y %H%M UT'),float(alts[alt_inx])))
            ax.set_title('\n'.join(txt),fontdict={'weight':'bold','size':'xx-large'})

            fname = '{0}_{1:03.0f}km_edens_map.png'.format(date.strftime('%Y%m%d_%H%MUT'),float(alts[alt_inx]))
            _filename = os.path.join(output_dir,fname)

            fig.savefig(_filename,bbox_inches='tight')
            plt.close()

            ################################################################################
            # Plot electron density profiles of endpoints for validation with IRI runs on CCMC
            # https://kauai.ccmc.gsfc.nasa.gov/instantrun/iri/
            if plot_profile_paths == 'all':
                for prof_key,profile in self.profiles.items():
                    fig     = plt.figure(figsize=(15,8))
                    nrows   = 2
                    ncols   = 1
                    ax_inx  = 0
                    pfxs    = ['tx','rx']
                    for pfx in pfxs:
                        ax_inx += 1
                        ax  = fig.add_subplot(nrows,ncols,ax_inx)

                        # Get latitude, altitude, and call of endpoint.
                        call     = profile.attrs['{!s}_call'.format(pfx)]
                        lat      = profile.attrs['{!s}_lat'.format(pfx)]
                        lon      = profile.attrs['{!s}_lon'.format(pfx)]
                        if lon > 180:
                            lon = lon - 360.

                        # Find closest lat/lon in currently calculated electron density array.
                        clst_lat_inx    = np.argmin(np.abs(lats-lat))
                        clst_lat        = lats[clst_lat_inx]

                        clst_lon_inx    = np.argmin(np.abs(lons-lon))
                        clst_lon        = lons[clst_lon_inx]

                        edp = edens[dInx,clst_lat_inx,clst_lon_inx,:]
                        ax.plot(alts,edp,marker='.')
                        ax.grid(True)

                        ax.set_xlabel('Altitude [km]')
                        ax.set_ylabel(r'IRI Electron Density [m$^{-3}$]')
                        actual      = '({:0.1f}\N{DEGREE SIGN} N, {:0.1f}\N{DEGREE SIGN} E)'.format(clst_lat,clst_lon)
                        requested   = '{!s} ({:0.1f}\N{DEGREE SIGN} N, {:0.1f}\N{DEGREE SIGN} E)'.format(call,lat,lon)

                        ax.set_title('Endpoint: {!s}\nActual Plotted: {!s}'.format(requested,actual))
                    
                    title   = []
                    title.append('IRI Endpoint Profiles')
                    title.append('{!s}'.format(date.strftime('%Y %b %d - %H%M UT')))
                    fig.text(0.5,1.00,'\n'.join(title),ha='center',va='bottom',fontdict={'weight':'bold','size':'large'})

                    fname = '{!s}_{!s}-{!s}_endPointProfiles'.format(date.strftime('%Y%m%d_%H%MUT'),
                            profile.attrs['tx_call'],profile.attrs['rx_call'])
                    _filename = os.path.join(output_dir,fname)

                    fig.tight_layout()
                    fig.savefig(_filename,bbox_inches='tight')
                    plt.close()

    def plot_maps_ortho(self,alt=250.,output_dir='output',figsize=(15,8),
            xlim=None,ylim=None,plot_profile_paths='all'):
        ds      = self.iri_dataset
        edens   = ds['electron_density'].values
        lats    = ds['lat'].values
        lons    = ds['lon'].values
        alts    = ds['alt'].values
        dates   = list(map(pd.to_datetime,ds['date'].values))

        alt_inx = np.argmin(np.abs(alts-alt))

        scale,ticks,cmap = calculate_scale(edens[:,:,:,alt_inx])
        bounds      = np.linspace(scale[0],scale[1],256)
        norm        = matplotlib.colors.BoundaryNorm(bounds,cmap.N)

        for dInx,date in tqdm.tqdm(enumerate(dates),desc='Plotting Maps',dynamic_ncols=True):
            tqdm.tqdm.write('Plotting: {!s}'.format(date))
            this_edens  = edens[dInx,:,:]

            fig     = plt.figure(figsize=figsize)
            ax      = fig.add_subplot(111,projection=ccrs.Orthographic(0,90))

    #        ax.coastlines(zorder=10,color='k')
    #        ax.add_feature(cartopy.feature.LAND)
    #        ax.add_feature(cartopy.feature.OCEAN)
    #        ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    #        ax.add_feature(cartopy.feature.RIVERS)
            ax.add_feature(cartopy.feature.COASTLINE)
            ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
            ax.set_title('')
            ax.gridlines(draw_labels=True)

            LONS, LATS  = np.meshgrid(lons,lats)
            data        = this_edens[:,:,alt_inx]
            pcoll       = ax.pcolormesh(LONS,LATS,data,cmap=cmap,norm=norm,transform=ccrs.PlateCarree())

            if plot_profile_paths == 'all':
                for prof_key,profile in self.profiles.items():
                    glons       = profile.attrs['glons']
                    tf          = glons > 180
                    glons[tf]   = glons[tf]-360.
                    glats       = profile.attrs['glats']
                    ax.plot(glons,glats,marker='o',color='k',lw=4,zorder=100,transform=ccrs.PlateCarree())

                    tx_call     = profile.attrs['tx_call']
                    tx_lat      = profile.attrs['tx_lat']
                    tx_lon      = profile.attrs['tx_lon']
                    if tx_lon > 180:
                        tx_lon = tx_lon - 360.

                    rx_call     = profile.attrs['rx_call']
                    rx_lat      = profile.attrs['rx_lat']
                    rx_lon      = profile.attrs['rx_lon']
                    if rx_lon > 180:
                        rx_lon = rx_lon - 360.

                    ax.scatter(tx_lon,tx_lat,marker='*',s=450,zorder=110,label=tx_call,ec='k',fc='red',transform=ccrs.PlateCarree())
                    ax.scatter(rx_lon,rx_lat,marker='*',s=450,zorder=110,label=rx_call,ec='k',fc='orange',transform=ccrs.PlateCarree())
                    fontdict = {'size':'x-large','weight':'bold'}
                    offset = 1.1
                    ax.text(tx_lon,tx_lat+offset,tx_call,fontdict=fontdict,ha='center',transform=ccrs.PlateCarree())
                    ax.text(rx_lon,rx_lat+offset,rx_call,fontdict=fontdict,ha='center',transform=ccrs.PlateCarree())

            # Plot Wave Source Locations
            if hasattr(self,'wave_list'):
                for wave in self.wave_list:
                    wlat    = wave['src_lat']
                    wlon    = wave['src_lon']
                    ax.scatter(wlon,wlat,marker='^',s=100,zorder=110,ec='k',fc='red',transform=ccrs.PlateCarree())


#            ax.set_xlim(xlim)
#            ax.set_ylim(ylim)

            ax.set_extent([-180, 180, 1, 90], ccrs.PlateCarree())
            cbar_label  = r'IRI Electron Density [m$^{-3}$]'
            cbar        = fig.colorbar(pcoll,orientation='vertical',shrink=0.65,pad=0.075,ticks=ticks)
            cbar.set_label(cbar_label,fontdict={'weight':'bold','size':'large'})

            # Plot Title
            txt = []
            txt.append('{0} - Alt: {1:.0f} km'.format(date.strftime('%d %b %Y %H%M UT'),float(alts[alt_inx])))
            ax.set_title('\n'.join(txt),fontdict={'weight':'bold','size':'xx-large'})

            fname = '{0}_{1:03.0f}km_edens_map.png'.format(date.strftime('%Y%m%d_%H%MUT'),float(alts[alt_inx]))
            _filename = os.path.join(output_dir,fname)

            fig.savefig(_filename,bbox_inches='tight')
            plt.close()

    def generate_wave(self,wave_list=None,advance_minutes=0.):
        """
        """
        ds      = self.iri_dataset
        lats    = ds['lat'].values
        lons    = ds['lon'].values
        alts    = ds['alt'].values
        dates   = list(map(pd.to_datetime,ds['date'].values))
        secs    = np.array([(x-dates[0]).total_seconds() for x in dates])

        SECS,LATS,LONS,ALTS  = np.meshgrid(secs,lats,lons,alts,indexing='ij')

        if wave_list is None:
            wave_list = []
            wave_list.append(dict(src_lat=60.,src_lon=-100.,amplitude=0.50,lambda_h=250,T_minutes=15))
#            wave_list.append(dict(src_lat=60.,src_lon= -80.,amplitude=0.50,lambda_h=250,T_minutes=15))

        for wave in wave_list:
            src_lat     = wave['src_lat']
            src_lon     = wave['src_lon']
            amplitude   = wave['amplitude']
            lambda_h    = wave['lambda_h']
            T_minutes   = wave['T_minutes']
            advance_minutes = wave.get('advance_minutes',0.)

            SRC_LATS    = np.ones_like(LATS) * src_lat
            SRC_LONS    = np.ones_like(LONS) * src_lon

            RANGES      = (Re + ALTS) * geopack.greatCircleDist(SRC_LATS,SRC_LONS,LATS,LONS)

            omega       = (2.*np.pi)/(T_minutes*60.)
            k_h         = (2.*np.pi)/lambda_h

            wave        = amplitude*np.cos(k_h*RANGES - omega*(SECS+advance_minutes*60.))
            ds          = ds + wave*ds

        ds.attrs['wave_list'] = str(wave_list)
        self.wave_list      = wave_list
        self.iri_dataset    = ds
