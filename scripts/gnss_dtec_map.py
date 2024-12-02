#!/usr/bin/env python
"""
This class will load and plot GNSS dTEC map.
"""
import os
import datetime
import tqdm
import h5py

import numpy as np
import pandas as pd
import xarray as xr
import scipy as sp
import scipy.stats

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import cartopy.crs as ccrs

mpl.rcParams['font.size']      = 12
mpl.rcParams['font.weight']    = 'bold'
mpl.rcParams['axes.grid']      = True
mpl.rcParams['grid.linestyle'] = ':'
mpl.rcParams['figure.figsize'] = np.array([15, 8])
mpl.rcParams['axes.xmargin']   = 0

# Region Dictionary
regions = {}
tmp     = {}
tmp['lon_lim']  = (-180.,180.)
tmp['lat_lim']  = ( -90., 90.)
regions['World']    = tmp

tmp     = {}
tmp['lon_lim']  = (-55.,-10)
tmp['lat_lim']  = (  40., 65.)
regions['Atlantic_Ocean']   = tmp

tmp     = {}
tmp['lon_lim']  = (-130.,-60.)
tmp['lat_lim']  = (  20., 55.)
regions['US']   = tmp

def get_bins(lim, bin_size):
    """ Helper function to split a limit into bins of the proper size """
    bins    = np.arange(lim[0], lim[1]+1*bin_size, bin_size)
    return bins

class GNSSdTECMap(object):
    def __init__(self,
            date      = datetime.datetime(2018,12,15,20),
            lat_lim   = (-90,90),
            lon_lim   = (-180,180),
            dlat      = 1.0,
            dlon      = 1.0,
            statistic = 'mean',
            data_dir  = None):

        if data_dir is None:
            data_dir    = os.path.join('data','gnss_dtec')

        self.data_dir  = data_dir
        self.date      = date
        self.lat_lim   = lat_lim
        self.lon_lim   = lon_lim
        self.dlat      = dlat
        self.dlon      = dlon
        self.statistic = statistic

        # Create a cache file name
        clst = []
        clst.append(date.strftime('%Y%m%d.%H%MUTC'))
        if lat_lim != (-90,90):
            clst.append('{!s}-{!s}lat'.format(*lat_lim))
        if lon_lim != (-180,180):
            clst.append('{!s}-{!s}lon'.format(*lon_lim))
        clst.append('{!s}dlat'.format(dlat))
        clst.append('{!s}dlon'.format(dlon))
        clst.append(statistic)
        clst.append('GNSSdTECmap.grid.nc')
        cfname  = '_'.join(clst)
        self.cfname = cfname
        cfpath  = os.path.join(data_dir,cfname)

        if os.path.exists(cfpath):
            self.gridded_dtec   = xr.load_dataarray(cfpath)
        else:
            self.load_data()
            self.grid_data()

    def load_data(self):
        data_dir        = self.data_dir
        date            = self.date

        date_str    = date.strftime('%Y%m%d')
        utc_str     = date.strftime('%H%M')
        #               20181215dTECmap2000.h5
        fname       = f'{date_str}dTECmap{utc_str}.h5'
        fpath       = os.path.join(data_dir,fname)
        if not os.path.exists(fpath):
            print('FILE NOT FOUND: {!s}'.format(fpath))
            return

        print('LOADING: {!s}'.format(fpath))
        self.fpath  = fpath
        with h5py.File(fpath,'r') as fl:
            data = {}
            data['lats'] = fl['la'][:]
            data['lons'] = fl['lo'][:]
            data['tid']  = fl['tid'][:]

        df      = pd.DataFrame(data)
        self.df = df

    def grid_data(self):
        df      = self.df

        lat_lim     = self.lat_lim
        lon_lim     = self.lon_lim
        dlat        = self.dlat
        dlon        = self.dlon
        statistic   = self.statistic

        lat_bins    = np.arange(lat_lim[0],lat_lim[1]+dlat,dlat)
        lon_bins    = np.arange(lon_lim[0],lon_lim[1]+dlon,dlon)
        lats        = lat_bins[:-1]
        lons        = lon_bins[:-1]

        tid_arr = np.zeros([len(lats),len(lons)])*np.nan
        xx      = df['lons'].values
        yy      = df['lats'].values
        vals    = df['tid'].values
        bins    = [lon_bins,lat_bins]

        print('GRIDDING dTEC DATA with scipy.stats.binned_statistic_2d()...')
        tic = datetime.datetime.now()
        res = sp.stats.binned_statistic_2d(xx,yy,vals,statistic=statistic,bins=bins)
        toc = datetime.datetime.now()
        print('   {!s}'.format(toc-tic))

        # Store in XArray
        xkey    = 'lon'
        ykey    = 'lat'

        crds    = {}
        crds[xkey]          = lons
        crds[ykey]          = lats

        attrs   = {}
        attrs['xkey']       = xkey
        attrs['ykey']       = ykey
        attrs['statistic']  = statistic
        bin_str = u'{:.1f}\N{DEGREE SIGN} lat x {:.1f}\N{DEGREE SIGN} lon bins'.format(dlat,dlon)
        attrs['label']      = '{!s} GNSS dTEC TECu\n'.format(statistic.title()) + bin_str
        
        arr                 = res.statistic
        gridded_dtec        = xr.DataArray(arr,crds,attrs=attrs,dims=[xkey,ykey])
        self.gridded_dtec   = gridded_dtec

        cfpath = os.path.join(self.data_dir,self.cfname)
        gridded_dtec.to_netcdf(cfpath)
        print('SAVING: {!s}'.format(cfpath))

    def plot_figure(self,png_fpath='output.png',figsize=(16,10),**kwargs):

        fig     = plt.figure(figsize=figsize)

        self.plot_map_ax(fig,subplot=(1,1,1))

        fig.tight_layout()
        fig.savefig(png_fpath,bbox_inches='tight')
        plt.close(fig)

    def plot_map_ax(self,fig,
            subplot             = (1,1,1),
            panel_rect          = None,
            plot_region         = 'World',
            title_size          = None,
            ticklabel_size      = None,
            label_size          = None,
            cbar_ticklabel_size = None,
            cbar_label_size     = None):

        if panel_rect is not None:
            ax = fig.add_axes(panel_rect,projection=ccrs.PlateCarree())
        else:
            ax = fig.add_subplot(*subplot,projection=ccrs.PlateCarree())

        # Plot Map #############################
        ax.coastlines(zorder=10,color='w')

        map_data        = self.gridded_dtec
        lon_key = map_data.attrs['xkey']
        lat_key = map_data.attrs['ykey']

        lons    = map_data[lon_key].values
        lats    = map_data[lat_key].values
        dlon    = lons[1] - lons[0]
        dlat    = lats[1] - lats[0]

        bin_str = '{:0.1f}'.format(dlat)+u'\N{DEGREE SIGN} lat x '+'{:0.1f}'.format(dlon)+u'\N{DEGREE SIGN} lon'
        map_data.name   = f'{self.statistic} dTEC TECu\n({bin_str})'

#        tf          = map_data < 1
#        map_n       = int(np.sum(map_data))
#        map_data    = np.log10(map_data)
#        map_data.values[tf] = 0
#        map_data.name   = 'log({})'.format(map_data.name)

        vmin    = -0.200
        vmax    =  0.200
        cntr    = map_data.plot.contourf(x=lon_key,y=lat_key,ax=ax,levels=30,cmap=mpl.cm.jet,vmin=vmin,vmax=vmax,cbar_kwargs=dict(shrink=0.65,pad=0.075))
#        cntr    = map_data.plot.pcolormesh(x=lon_key,y=lat_key,ax=ax,cmap=mpl.cm.jet,vmin=vmin,vmax=vmax)
        cax     = cntr.colorbar.ax
        cax.set_ylabel(map_data.name)
        if cbar_ticklabel_size is not None:
            for ytl in cax.get_yticklabels():
                ytl.set_size(cbar_ticklabel_size)
        if cbar_label_size is not None:
            cax.set_ylabel(map_data.name,fontdict={'size':cbar_label_size})

        lweight = mpl.rcParams['axes.labelweight']
        lsize   = mpl.rcParams['axes.labelsize']
        fdict   = {'weight':lweight,'size':lsize}
        if label_size is not None:
            fdict.update({'size':label_size})

        if plot_region is not None:
            rgn         = regions.get(plot_region)
            lat_lim     = rgn.get('lat_lim')
            lon_lim     = rgn.get('lon_lim')

            ax.set_xlim(lon_lim)
            ax.set_ylim(lat_lim)

        if ticklabel_size is not None:
            for ttl in ax.get_xticklabels():
                ttl.set_size(ticklabel_size)

            for ttl in ax.get_yticklabels():
                ttl.set_size(ticklabel_size)

        date_str = self.date.strftime('%Y %b %d')
        title_fd = {'weight':'bold'}
        if title_size is not None:
            title_fd.update({'size':title_size})
        ax.set_title(date_str,fontdict=title_fd)
        return ax

if __name__ == '__main__':
    output_dir = os.path.join('output','GNSSdTECmap')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    date    = datetime.datetime(2018,12,15,20)
    hsp     = GNSSdTECMap(date)

    png_fname   = date.strftime('%Y%m%d.%H%M')+'_GNSSdTECmap.png'
    png_fpath   = os.path.join(output_dir,png_fname)
    hsp.plot_figure(png_fpath=png_fpath)
    print(png_fpath)
