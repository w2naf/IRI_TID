#!/usr/bin/env python
import os
import matplotlib

from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

from ionolib import geopack

# Geographiclib - https://geographiclib.sourceforge.io/Python/2.0/
from geographiclib.geodesic import Geodesic
geod = Geodesic.WGS84

Re = 6371.

def plot_maps(profiles={},output_dir='output',fname='geopack_test.png',
        figsize=(15,8),xlim=None,ylim=None,plot_profile_paths='all'):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fig     = plt.figure(figsize=figsize)
    ax      = fig.add_subplot(111,projection=ccrs.PlateCarree())

#    ax.coastlines(zorder=10,color='k')
#    ax.add_feature(cartopy.feature.LAND)
#    ax.add_feature(cartopy.feature.OCEAN)
#    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
#    ax.add_feature(cartopy.feature.RIVERS)

    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
    ax.set_title('')
    ax.gridlines(draw_labels=True)

    for prof_key,profile in profiles.items():
        glons       = profile['glons']
        tf          = glons > 180
        glons[tf]   = glons[tf]-360.
        glats       = profile['glats']
        line        = ax.plot(glons,glats,marker='o',ls=' ',lw=2,zorder=100,label=prof_key)

        call_0     = profile['call_0']
        lat_0      = profile['lat_0']
        lon_0      = profile['lon_0']
        if lon_0 > 180:
            lon_0 = lon_0 - 360.

        call_1     = profile['call_1']
        lat_1      = profile['lat_1']
        lon_1      = profile['lon_1']
        if lon_1 > 180:
            lon_1 = lon_1 - 360.

        ax.scatter(lon_0,lat_0,marker='*',s=450,zorder=110,ec='k',fc='red')
        ax.scatter(lon_1,lat_1,marker='*',s=450,zorder=110,ec='k',fc='orange')
        fontdict = {'size':'x-large','weight':'bold'}
        offset = 1.1
        ax.text(lon_0,lat_0+offset,call_0,fontdict=fontdict,ha='center')
        ax.text(lon_1,lat_1+offset,call_1,fontdict=fontdict,ha='center')

    ax.legend(loc='lower right')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    _filename = os.path.join(output_dir,fname)

    fig.savefig(_filename,bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    call_0      = 'W2NAF'
    lat_0       =   41.335116
    lon_0       =  -75.600692

    call_1      = 'WWV'
    lat_1       =   40.6683
    lon_1       = -105.0384

    range_step  = 100.

    methods     = []
    methods.append('geopack')
    methods.append('geographiclib')

    profiles = {}
    for method in methods:
        prfl = profiles[method] = {}

        if method == 'geopack':
            # Determine the ranges and azimuth along the profile path.
            dist        = Re * geopack.greatCircleDist(lat_0,lon_0,lat_1,lon_1)
            ranges      = np.arange(0,dist,range_step)
            az          = geopack.greatCircleAzm(lat_0,lon_0,lat_1,lon_1)

            # Calculate the lats/lons to be used along the profile path (field points).
            flats_flons   = []
            for x in ranges:
                tmp = geopack.calcDistPnt(lat_0,lon_0,origAlt=0,az=az,dist=x,el=0)
                flats_flons.append([tmp['distLat'],tmp['distLon']])
            flats_flons         = np.array(flats_flons)
            glats   = flats_flons[:,0]
            glons   = flats_flons[:,1]

        elif method == 'geographiclib':
            invl    = geod.InverseLine(lat_0,lon_0,lat_1,lon_1)
            dist    = invl.s13*1e-3   # Distance in km
            az      = invl.azi1

            ranges  = np.arange(0,dist,range_step)

            glats   = []
            glons   = []
            for x in ranges:
                s   = min(x*1e3,invl.s13) # invl.s13 is the total line distance in m
                tmp = invl.Position(s,Geodesic.STANDARD)
                invl_dist   = tmp['s12']
                glat        = tmp['lat2']
                glon        = tmp['lon2']
                azi2        = tmp['azi2']

                glats.append(glat)
                glons.append(glon)
            
        prfl['call_0']  = call_0
        prfl['lat_0']   = lat_0
        prfl['lon_0']   = lon_0

        prfl['call_1']  = call_1
        prfl['lat_1']   = lat_1
        prfl['lon_1']   = lon_1
        
        prfl['dist']    = dist
        prfl['az']      = az
        prfl['ranges']  = ranges
        prfl['glats']   = np.array(glats)
        prfl['glons']   = np.array(glons)

## World
#xlim    = (-180,180)
#ylim    = (-90,90)

## CONUS
#xlim    = (-130,-56)
#ylim    = (20,55)

# CONUS + Canada
prmd    = {}
prmd['xlim']        = (-130,-56)
prmd['ylim']        =  (20, 80)
prmd['profiles']    = profiles
plot_maps(**prmd)

import ipdb; ipdb.set_trace()
