#!/usr/bin/env python
import os
import shutil
import datetime
import pickle

import matplotlib

from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon

from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import scipy.interpolate
import pandas as pd

import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import pyiri2016

try:
    from . import geopack
except:
    import geopack


Re = 6371 # Earth Radius in km

def calculate_scale(data,stddevs=2.,lim=(-125000,125000)):
    if lim == 'auto':
        mean    = np.nanmean(np.abs(data))
        std     = np.nanstd(np.abs(data))
        lim     = mean + stddevs*std
        scale   = (-lim,lim)
    else:
        scale   = lim

    ticks   = None
    cmap    = matplotlib.cm.jet

    return scale,ticks,cmap

class iono_3d(object):
    def __init__(self,date,
            hgt_0 =   60., hgt_1 =  560., hgt_step=1.0,
            lat_0 =   32., lat_1 =   80., lat_step=1.00,
            lon_0 = -100., lon_1 =  -40., lon_step=1.00):

        dhour       = date.hour + date.minute/60.

        lats    = np.arange(lat_0,lat_1,lat_step)
        lons    = np.arange(lon_0,lon_1,lon_step)
        alts    = np.arange(hgt_0,hgt_1,hgt_step)

        self.lats       = lats
        self.lons       = lons % 360. # Adjust so lons are between 0 and 360 deg.
        self.alts       = alts
        self.lat_step   = lat_step
        self.lon_step   = lon_step
        self.alt_step   = hgt_step
        self.date       = date
        self.iri_date   = date  # The IRI background date/time may be different from the
                                # one used by an artificial wave.
        self.profiles   = {}

    def generate_tx_rx_profile(self,tx_lat,tx_lon,rx_lat,rx_lon,range_step=50.,tx_call='None',rx_call='None',date=None):
        dist    = Re * geopack.greatCircleDist(tx_lat,tx_lon,rx_lat,rx_lon)
        ranges  = np.arange(0,dist,range_step)

        if date is None: date = self.date
        dhour   = date.hour + date.minute/60.

        az          = geopack.greatCircleAzm(tx_lat,tx_lon,rx_lat,rx_lon)
        flat_flon   = []
        for x in ranges:
            tmp = geopack.calcDistPnt(tx_lat,tx_lon,origAlt=0,az=az,dist=x,el=0)
            flat_flon.append([tmp['distLat'],tmp['distLon']])

        flat_flon       = np.array(flat_flon)

        shape           = flat_flon.shape
        flat_flon       = flat_flon.flatten()
        tf              = flat_flon < 0
        flat_flon[tf]   = 360. + flat_flon[tf]
        flat_flon.shape = shape

        dict_key = '{}-{}'.format(tx_call,rx_call)

        e_density   = np.zeros([len(ranges),len(self.alts)])
        for inx,(lat,lon) in enumerate(flat_flon):
            print('IRI: {!s} Lat: {:0.2f} Lon: {:0.2f}'.format(date,lat,lon))

            #  jmag    --> 0: geographic; 1: geomagnetic
            #  iut     --> 0: for LT;     1: for UT
            #  htecmax --> 0: no TEC; otherwise: upper boundary for integral [km]
            result = pyiri2016.IRI2016Profile(altlim=[self.alts.min(),self.alts.max()],
                    altstp=self.alt_step,lat=lat,lon=lon,
                    year=date.year,month=date.month,dom=date.day,hour=dhour,
                    jmag=0,iut=0,htecmax=0,verbose=False)

            hei_df          = result.heiProfile
            e_density[inx,:]    = hei_df['Ne']

        # Calculate range from start as an angle [radians]
        # Computed the same way as in raydarn fortran code.
        field_lats  = flat_flon[:,0]
        field_lons  = flat_flon[:,1]
        edensTHT    = np.arccos( np.cos(field_lats[0]*np.pi/180.)*np.cos(field_lats*np.pi/180.)* \
                            np.cos((field_lons - field_lons[0])*np.pi/180.) \
                    + np.sin(field_lats[0]*np.pi/180.)*np.sin(field_lats*np.pi/180.))

        # Set initial edensTHT = 0
        edensTHT[0] = 0

        profl = {}
        profl['e_density']  = e_density
        profl['glats']      = flat_flon[:,0]
        profl['glons']      = flat_flon[:,1]
        profl['alts']       = self.alts
        profl['ranges']     = ranges
        profl['dip']        = e_density*0.
        profl['edensTHT']   = edensTHT

        # Base filename to be used with this profile.
        date_code   = self.date.strftime('%Y%m%d_%H%MUT')
        fname_base  = '{date_code}_{tx_call}_{rx_call}'.format(
                date_code=date_code,tx_call=tx_call,rx_call=rx_call)

        # Add extra meta-data.
        profl['date']       = self.date
        profl['tx_call']    = tx_call
        profl['rx_call']    = rx_call
        profl['fname_base'] = fname_base

        self.profiles[dict_key] = profl
        return self.profiles

    def run_iri(self):
        for lat_inx,lat in enumerate(lats):
            for lon_inx,lon in enumerate(lons):

                print('IRI: {!s} Lat: {:0.2f} Lon: {:0.2f}'.format(date,lat,lon))

                #  jmag    --> 0: geographic; 1: geomagnetic
                #  iut     --> 0: for LT;     1: for UT
                #  htecmax --> 0: no TEC; otherwise: upper boundary for integral [km]
                import ipdb; ipdb.set_trace()
                result = pyiri2016.IRI2016Profile(altlim=[hgt_0,hgt_1], altstp=hgt_step, 
                    year=date.year,month=date.month,dom=date.day,hour=dhour,
                    lat=lat, lon=lon, jmag=0,iut=0,htecmax=0,verbose=False)

                hei_df = result.heiProfile

                edens[lat_inx,lon_inx,:]    = hei_df['Ne']
#                dip[lat_inx,lon_inx,0]      = oarr[24] # Dip [deg]
#                dip[lat_inx,lon_inx,1]      = oarr[26] # MODIFIED DIP LATITUDE

    def plot_profiles(self,keys=None,output_dir='output',filename=None,figsize=(10,6)):

        if keys is None:
            keys = self.profiles.keys()
        
        paths = []
        for key in keys:
            profl   = self.profiles[key]
            edens   = profl['e_density']
            lats    = profl['glats']
            lons    = profl['glons']
            alts    = profl['alts']
            ranges  = profl['ranges']

            range_step  = np.mean(np.diff(ranges))
            alt_step    = np.mean(np.diff(alts))

            fig  = plt.figure(figsize=figsize)
            ax   = fig.add_subplot(111)

            scale,ticks,cmap = calculate_scale(edens)
            bounds  = np.linspace(scale[0],scale[1],256)
            norm    = matplotlib.colors.BoundaryNorm(bounds,cmap.N)
            
            pcoll = ax.pcolorfast(ranges,alts,edens.T,cmap=cmap,norm=norm)

            ax.set_xlim(ranges.min(), ranges.max())
            ax.set_ylim(alts.min(), alts.max())

            ax.set_xlabel('Range [km]')
            ax.set_ylabel('Altitude [km]')

            cbar    = fig.colorbar(pcoll,orientation='vertical',shrink=0.60,pad=.10,ticks=ticks)
            txt     = r'IRI Electron Density [m$^{-3}$]'
            cbar.set_label(txt)

            txt = []
            txt.append('IRI Electron Density')
            txt.append('{0} {1}'.format(key,self.date.strftime('%d %b %Y %H%M UT')))
            ax.set_title('\n'.join(txt))

            fig.tight_layout()

            if filename is None:
                _filename = os.path.join(output_dir,'{0}_profile.png'.format(profl['fname_base']))
            else:
                _filename = os.path.join(output_dir,filename)

            fig.savefig(_filename,bbox_inches='tight')
            plt.close()
            paths.append(_filename)
        return paths

    def plot_maps(self,keys=None,alt=250.,output_dir='output',filename=None,figsize=(10,8)):

        if keys is None:
            keys = self.profiles.keys()
        
        for key in keys:
            profl   = self.profiles[key]
            edens   = profl['e_density']
            lats    = profl['glats']
            lons    = profl['glons']
            alts    = profl['alts']
            ranges  = profl['ranges']

            fig     = plt.figure(figsize=figsize)
            ax      = fig.add_subplot(111,projection=ccrs.PlateCarree())
            ax_only = False

    #        ax.coastlines(zorder=10,color='k')
    #        ax.add_feature(cartopy.feature.LAND)
    #        ax.add_feature(cartopy.feature.OCEAN)
            ax.add_feature(cartopy.feature.COASTLINE)
            ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
    #        ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    #        ax.add_feature(cartopy.feature.RIVERS)
            ax.set_title('')


            alt_inx = np.argmin(np.abs(self.alts-alt))

            try:
                scale,ticks,cmap = calculate_scale(edens[:,:,alt_inx])
            except:
                import ipdb; ipdb.set_trace()
            bounds      = np.linspace(scale[0],scale[1],256)
            norm        = matplotlib.colors.BoundaryNorm(bounds,cmap.N)

            LONS, LATS  = np.meshgrid(self.lons,self.lats)
            data        = self.edens[:,:,alt_inx]
            pcoll       = ax.pcolormesh(LONS,LATS,data,cmap=cmap,norm=norm)

            cbar_label  = r'IRI Electron Density [m$^{-3}$]'
            cbar    = fig.colorbar(pcoll,orientation='vertical',shrink=0.75,pad=0.075,ticks=ticks)
            cbar.set_label(cbar_label,fontdict={'weight':'bold','size':'large'})

            if plot_title:
                txt = []
        #        txt.append('IRI Electron Density')
                txt.append('{0} - Alt: {1:.0f} km'.format(self.date.strftime('%d %b %Y %H%M UT'),float(self.alts[alt_inx])))
                ax.set_title('\n'.join(txt),fontdict={'weight':'bold','size':'large'})

            fig.tight_layout()

            if filename is None:
                fname = '{0}_{1:03.0f}km_edens_map.png'.format(self.date.strftime('%Y%m%d_%H%MUT'),float(self.alts[alt_inx]))
                _filename = os.path.join(output_dir,fname)
            else:
                _filename = os.path.join(output_dir,filename)

            fig.savefig(_filename,bbox_inches='tight')
            plt.close()

    def generate_wave(self,wave_list=None,sTime=None,currTime=None,keys=None):
        """
        sTime: Set to self.iri_date if None.
        """
        if wave_list is None:
            wave_list = []
            wave_list.append(dict(src_lat=60.,src_lon=-100.,amplitude=0.50,lambda_h=250,T_minutes=15))
            wave_list.append(dict(src_lat=60.,src_lon=-80.,amplitude=0.50,lambda_h=250,T_minutes=15))

        if keys is None:
            keys = self.profiles.keys()
        
        paths = []
        for key in keys:
            profl   = self.profiles[key]
            edens   = profl['e_density']
            lats    = profl['glats']
            lons    = profl['glons']
            alts    = profl['alts']
            ranges  = profl['ranges']

            range_step   = np.mean(np.diff(ranges))
            alt_step     = np.mean(np.diff(alts))

            rng_from_src = np.zeros(edens.shape)

            edens_0 = edens.copy()
            for wave in wave_list:
                src_lat     = wave['src_lat']
                src_lon     = wave['src_lon']
                amplitude   = wave['amplitude']
                lambda_h    = wave['lambda_h']
                T_minutes   = wave['T_minutes']

                src_lats = np.ones_like(lats) * src_lat
                src_lons = np.ones_like(lons) * src_lon

                # rng_from_src = (Re + alts) * geopack.greatCircleDist(src_lats,src_lons,lats,lons)
                realts      = np.broadcast_to((Re + alts),edens.shape)
                dists       = geopack.greatCircleDist(src_lats,src_lons,lats,lons)
                dists       = np.broadcast_to(dists,edens.shape[::-1]).T
                rng_from_src = realts * dists

                if sTime is None: sTime = self.iri_date

                if (T_minutes is None) or (currTime is None):
                    omega = 0.
                    t = 0.
                else:
                    omega = (2.*np.pi)/(T_minutes*60.)
                    t = (currTime - sTime).total_seconds()
                    self.date = currTime

                k_h         = 2.*np.pi/lambda_h
                wave        = edens_0 * amplitude*np.sin(k_h*rng_from_src - omega*t)
                edens      += wave
            profl['e_density'] = edens

        self.wave_list = wave_list
