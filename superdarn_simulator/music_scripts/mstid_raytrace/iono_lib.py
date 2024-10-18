#!/usr/bin/env python
import os
import shutil
import datetime
import pickle

import multiprocessing

import matplotlib
matplotlib.use('Agg')

from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd

from davitpy.models import *
from davitpy import utils
from davitpy import pydarn
from naf_tools import prepare_output_dirs

def calculate_scale(data,stddevs=2.,lim='auto'):
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

def make_dirs(base='data/event',clear_output_dirs=False,make_dirs=True):
    dirs_list = []
    dirs_list.append('000_general_output')
    dirs_list.append('background_edens_runfiles')
    dirs_list.append('background_edens_raw')
    dirs_list.append('background_edens_perturbed')
    dirs_list.append('edens_profl_dats')
    dirs_list.append('edens_alt_p')
    dirs_list.append('edens_profl_figs')
    dirs_list.append('edens_maps')
    dirs_list.append('rt_data')
    dirs_list.append('rt_profl_figs')
    dirs_list.append('rt_maps')

    if clear_output_dirs:
        try:
            if os.path.isdir(base) and not os.path.islink(base):
                shutil.rmtree(base)
            elif os.path.islink(base):
                os.remove(base)
        except:
            pass

    dirs_dict = {}
    for dr in dirs_list:
        dirs_dict[dr] = os.path.join(base,dr)

    if make_dirs:
        prepare_output_dirs(dirs_dict,clear_output_dirs=clear_output_dirs)

    dirs_dict['base'] = base
    return dirs_dict

#            hgt_0 =   60., hgt_1 =  560., hgt_step=1.0,
#            lat_0 =   20., lat_1 =   60., lat_step=0.1,
#            lon_0 = -110., lon_1 =  -60., lon_step=0.1,
class iono_3d(object):
    def __init__(self,date,
            hgt_0 =   60., hgt_1 =  560., hgt_step=1.0,
            lat_0 =   32., lat_1 =   80., lat_step=1.00,
            lon_0 = -100., lon_1 =  -40., lon_step=1.00,
            load_from_pickle=False):
#            hgt_0 =   60., hgt_1 =  560., hgt_step=1.0,
#            lat_0 =   30., lat_1 =   60., lat_step=0.1,
#            lon_0 = -110., lon_1 =  -70., lon_step=0.1,
#            load_from_pickle=False):

        # Inputs
        jf = [True]*50
        jf[2:6] = [False]*4
        jf[20] = False
        jf[22] = False
        jf[27:30] = [False]*3
        jf[32] = False
        jf[34] = False
        jmag = 0.
        iyyyy = date.year
        mmdd = 100*date.month+date.day
        dhour = date.hour + date.minute/60. + 25

        lats    = np.arange(lat_0,lat_1,lat_step)
        lons    = np.arange(lon_0,lon_1,lon_step)
        alts    = np.arange(hgt_0,hgt_1,hgt_step)


        edens    = np.zeros((lats.size,lons.size,alts.size),dtype=np.float32)
        dip      = np.zeros((lats.size,lons.size,2),dtype=np.float32)
        for lat_inx,lat in enumerate(lats):
            for lon_inx,lon in enumerate(lons):

                print 'IRI: {!s} Lat: {:0.2f} Lon: {:0.2f}'.format(date,lat,lon)

                oarr = np.zeros(100)
                # Call fortran subroutine
                outf,oarr = iri.iri_sub(jf,jmag,lat,lon,iyyyy,mmdd,dhour,hgt_0,hgt_1,hgt_step,oarr)

                edens[lat_inx,lon_inx,:]    = outf[0,0:alts.size]
                dip[lat_inx,lon_inx,0]      = oarr[24] # Dip [deg]
                dip[lat_inx,lon_inx,1]      = oarr[26] # MODIFIED DIP LATITUDE

        self.lats       = lats
        self.lons       = lons % 360. # Adjust so lons are between 0 and 360 deg.
        self.alts       = alts
        self.lat_step   = lat_step
        self.lon_step   = lon_step
        self.alt_step   = hgt_step
        self.edens      = edens
        self.date       = date
        self.iri_date   = date  # The IRI background date/time may be different from the
                                # one used by an artificial wave.
        self.dip        = dip
        self.profiles   = {}

        #*  **JF switches to turn off/on (True/False) several options**
        #
        #*  [0] :    True
        #    * Ne computed            
        #    * Ne not computed
        #*  [1] :    True
        #    * Te, Ti computed        
        #    * Te, Ti not computed
        #*  [2] :    True
        #    * Ne & Ni computed       
        #    * Ni not computed
        #*  [3] :    False
        #    * B0 - Table option      
        #    * B0 - other models jf[30]
        #*  [4] :    False
        #    * foF2 - CCIR            
        #    * foF2 - URSI
        #*  [5] :    False
        #    * Ni - DS-95 & DY-85     
        #    * Ni - RBV-10 & TTS-03
        #*  [6] :    True
        #    * Ne - Tops: f10.7<188   
        #    * f10.7 unlimited          
        #*  [7] :    True
        #    * foF2 from model        
        #    * foF2 or NmF2 - user input
        #*  [8] :    True
        #    * hmF2 from model        
        #    * hmF2 or M3000F2 - user input
        #*  [9] :    True
        #    * Te - Standard          
        #    * Te - Using Te/Ne correlation
        #* [10] :    True
        #    * Ne - Standard Profile  
        #    * Ne - Lay-function formalism
        #* [11] :    True
        #    * Messages to unit 6     
        #    * to meesages.text on unit 11
        #* [12] :    True
        #    * foF1 from model        
        #    * foF1 or NmF1 - user input
        #* [13] :    True
        #    * hmF1 from model        
        #    * hmF1 - user input (only Lay version)
        #* [14] :    True
        #    * foE  from model        
        #    * foE or NmE - user input
        #* [15] :    True
        #    * hmE  from model        
        #    * hmE - user input
        #* [16] :    True
        #    * Rz12 from file         
        #    * Rz12 - user input
        #* [17] :    True
        #    * IGRF dip, magbr, modip 
        #    * old FIELDG using POGO68/10 for 1973
        #* [18] :    True
        #    * F1 probability model   
        #    * critical solar zenith angle (old)
        #* [19] :    True
        #    * standard F1            
        #    * standard F1 plus L condition
        #* [20] :    False
        #    * ion drift computed     
        #    * ion drift not computed
        #* [21] :    True
        #    * ion densities in %     
        #    * ion densities in m-3
        #* [22] :    False
        #    * Te_tops (Aeros,ISIS)   
        #    * Te_topside (TBT-2011)
        #* [23] :    True
        #    * D-region: IRI-95       
        #    * Special: 3 D-region models
        #* [24] :    True
        #    * F107D from APF107.DAT  
        #    * F107D user input (oarr[41])
        #* [25] :    True
        #    * foF2 storm model       
        #    * no storm updating
        #* [26] :    True
        #    * IG12 from file         
        #    * IG12 - user
        #* [27] :    False
        #    * spread-F probability   
        #    * not computed
        #* [28] :    False
        #    * IRI01-topside          
        #    * new options as def. by JF[30]
        #* [29] :    False
        #    * IRI01-topside corr.    
        #    * NeQuick topside model
        #* [28,29]:
        #    * [t,t] IRIold, 
        #    * [f,t] IRIcor, 
        #    * [f,f] NeQuick, 
        #    * [t,f] Gulyaeva
        #* [30] :    True
        #    * B0,B1 ABT-2009     
        #    * B0 Gulyaeva h0.5
        #* [31] :    True
        #    * F10.7_81 from file     
        #    * PF10.7_81 - user input (oarr[45])
        #* [32] :    False
        #    * Auroral boundary model on
        #    * Auroral boundary model off
        #* [33] :    True
        #    * Messages on            
        #    * Messages off
        #* [34] :    False
        #    * foE storm model        
        #    * no foE storm updating
        #* [..] :    ....
        #* [50] :    ....

#        self.outf   = outf
#        self.oarr   = oarr

    def plot_profiles(self,keys=None,output_dir='data',filename=None,figsize=(10,6)):

        if keys is None:
            keys = self.profiles.keys()
        
        keys = np.array(keys)
        
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

    def plot_map(self,figsize=(10,8),alt=250.,radar=None,fov=None,
            ax = None,
            plot_title  = True,
            output_dir='data',filename=None,
            basemap_dict =  {'lat_0':60., 'lon_0':-105, 'width':6500e3, 'height':6500e3, 'lat_ts':50., 'projection':'stere','resolution':'l'}):

        if ax is None:
            fig     = plt.figure(figsize=figsize)
            ax      = fig.add_subplot(111)
            ax_only = False
        else:
            fig     = ax.get_figure()
            ax_only = True


        ax_info = {}

#        m = Basemap(ax=ax,**basemap_dict)
        m = utils.plotUtils.mapObj(ax=ax,datetime=self.date,fillContinents='None',
                fillLakes='None',**basemap_dict)

        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.),color='k',labels=[False,True,False,False])
        m.drawmeridians(np.arange(-180.,181.,20.),color='k',labels=[True,False,False,True])
#        m.drawparallels(np.arange(-80.,81.,20.),color='k',labels=[False,False,False,False])
#        m.drawmeridians(np.arange(-180.,181.,20.),color='k',labels=[False,False,False,False])
        m.drawcoastlines(color='1.',linewidth=4)
        m.drawcoastlines(linewidth=1)
        m.drawmapboundary(fill_color='w')

        alt_inx = np.argmin(np.abs(self.alts-alt))

        scale,ticks,cmap = calculate_scale(self.edens[:,:,alt_inx])
        bounds      = np.linspace(scale[0],scale[1],256)
        norm        = matplotlib.colors.BoundaryNorm(bounds,cmap.N)

        LONS, LATS  = np.meshgrid(self.lons,self.lats)
        data        = self.edens[:,:,alt_inx]
        pcoll       = m.pcolormesh(LONS,LATS,data,cmap=cmap,norm=norm,latlon=True)

        cbar_label  = r'IRI Electron Density [m$^{-3}$]'
        ax_info['cbar_pcoll']   = pcoll
        ax_info['cbar_ticks']   = ticks 
        ax_info['cbar_label']   = cbar_label 

        if not ax_only:
            cbar    = fig.colorbar(pcoll,orientation='vertical',shrink=0.75,pad=0.075,ticks=ticks)
            cbar.set_label(cbar_label,fontdict={'weight':'bold','size':'large'})

        #Overlay radar field of view.
        if radar is not None:
            if fov is None:
                site    = pydarn.radar.radStruct.site(code=radar,dt=self.date)
                fov     = pydarn.radar.radFov.fov(site=site,model='IS',coords='geo')

            verts = []
            data  = []
            for beam_inx,beam in enumerate(fov.beams):
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
                    
            pcoll = PolyCollection(np.array(verts),facecolor=(1,1,1,0),edgecolors='black',linewidths=None,closed=False,cmap=cmap,norm=norm,zorder=200)
            ax.add_collection(pcoll,autolim=False)

            pydarn.plotting.overlayRadar(m,codes=[radar],markerSize=30,fontSize=20,yOffset=-15)

        if plot_title:
            txt = []
    #        txt.append('IRI Electron Density')
            txt.append('{0} - Alt: {1:.0f} km'.format(self.date.strftime('%d %b %Y %H%M UT'),float(self.alts[alt_inx])))
            ax.set_title('\n'.join(txt),fontdict={'weight':'bold','size':'large'})

        if not ax_only:
            fig.tight_layout()

            if filename is None:
                fname = '{0}_{1:03.0f}km_edens_map.png'.format(self.date.strftime('%Y%m%d_%H%MUT'),float(self.alts[alt_inx]))
                _filename = os.path.join(output_dir,fname)
            else:
                _filename = os.path.join(output_dir,filename)

            fig.savefig(_filename,bbox_inches='tight')
            plt.close()

        return ax_info

    def plot_lon(self,lon=-95,figsize=(10,8)):

        fig  = plt.figure(figsize=figsize)
        ax   = fig.add_subplot(111)


        # Figure out which longitude in the array is closest to the one of interest.
        if lon < 0:
            tmp_lon = lon + 360.
        else:
            tmp_lon = lon

        # Make handle the fact that longitude is cyclical.
        tmp_lons = np.array(self.lons)
        tmp_lons[tmp_lons < 0] = tmp_lons[tmp_lons < 0] + 360.
        lon_inx = np.argmin( np.abs(tmp_lons-tmp_lon)%360 )

        # Get values of lat/alt_step
        lat_step = self.lat_step
        alt_step = self.alt_step

        verts = []
        data  = []
        for lat_inx,lat in enumerate(self.lats):
            for alt_inx,alt in enumerate(self.alts):
                x1,y1 = lat-lat_step/2.,alt-alt_step/2.
                x2,y2 = lat+lat_step/2.,alt-alt_step/2.
                x3,y3 = lat+lat_step/2.,alt+alt_step/2.
                x4,y4 = lat-lat_step/2.,alt+alt_step/2.
                verts.append(((x1,y1),(x2,y2),(x3,y3),(x4,y4),(x1,y1)))
                data.append(self.edens[lat_inx,lon_inx,alt_inx])

        data = np.array(data)
        scale,ticks,cmap = self.__calculate_scale(data)
        bounds  = np.linspace(scale[0],scale[1],256)
        norm    = matplotlib.colors.BoundaryNorm(bounds,cmap.N)

        pcoll = PolyCollection(np.array(verts),edgecolors='face',linewidths=0,closed=False,cmap=cmap,norm=norm)
        pcoll.set_array(data)
        ax.add_collection(pcoll,autolim=False)

        ax.set_xlim(self.lats.min(), self.lats.max())
        ax.set_ylim(self.alts.min(), self.alts.max())

        ax.set_xlabel('Latitude [deg]')
        ax.set_ylabel('Altitude [km]')

        cbar    = fig.colorbar(pcoll,orientation='vertical',shrink=0.60,pad=.10,ticks=ticks)
        txt     = r'IRI Electron Density [m$^{-3}$]'
        cbar.set_label(txt)

        txt = []
        txt.append('IRI Electron Density')
        txt.append('{0} - GLon: {1:.1f} deg'.format(self.date.strftime('%d %b %Y %H%M UT'),float(self.lons[lon_inx])))
        ax.set_title('\n'.join(txt))

        fig.tight_layout()
        fig.savefig(os.path.join('data','edens_long.png'),bbox_inches='tight')
        plt.close()

    def generate_profile(self, field_lats, field_lons, ranges, dict_key,debug=False):
        # Electron Density Interpolation ###############################################
        field_lats = np.array(field_lats)
        field_lons = np.array(field_lons)
        field_alts = np.array(self.alts)
        if field_lats.shape != field_lons.shape:
            print 'ERROR: field_lats and field_lons are not the same shape.'
            return None

        # Generate mesh grid for field points.
        field_shape = (len(field_lats),field_alts.size)
        f_lat_mesh = np.zeros(field_shape,dtype=np.float32)
        for inx in range(field_shape[1]):
            f_lat_mesh[:,inx] = field_lats

        f_lon_mesh = np.zeros(field_shape,dtype=np.float32)
        for inx in range(field_shape[1]):
            f_lon_mesh[:,inx] = field_lons

        f_alt_mesh = np.zeros(field_shape,dtype=np.float32)
        for inx in range(field_shape[0]):
            f_alt_mesh[inx,:] = field_alts

        points  = (self.lats, self.lons, self.alts)
        values  = self.edens
        xi      = np.array((f_lat_mesh.flatten(),f_lon_mesh.flatten(),f_alt_mesh.flatten())).T

        if debug:
            t0 = datetime.datetime.now()
            print 'Starting interpolation: ', t0
        field_profile = sp.interpolate.interpn(points, values, xi, method='nearest',bounds_error=False,fill_value=None)
        if debug:
            t1 = datetime.datetime.now()
            print 'Done interpolating: ', t1
            print 'Total time: ',(t1 - t0)

        field_profile.shape = field_shape
        # End Electron Density Interpolation ###########################################

        # Calculate range from start as an angle [radians]
        # Computed the same way as in raydarn fortran code.
        edensTHT    = np.arccos( np.cos(field_lats[0]*np.pi/180.)*np.cos(field_lats*np.pi/180.)* \
                            np.cos((field_lons - field_lons[0])*np.pi/180.) \
                    + np.sin(field_lats[0]*np.pi/180.)*np.sin(field_lats*np.pi/180.))

        # Dip Parameter Interpolation ################################################## 
        points  = (self.lats ,self.lons)
        xi      = (field_lats, field_lons)

        values  = self.dip[:,:,0]
        dip_0   = sp.interpolate.interpn(points, values, xi, method='nearest',bounds_error=False,fill_value=None)

        values  = self.dip[:,:,1]
        dip_1   = sp.interpolate.interpn(points, values, xi, method='nearest',bounds_error=False,fill_value=None)

        profl_dip = np.zeros([field_lats.size,2],dtype=np.float32)
        profl_dip[:,0] = dip_0
        profl_dip[:,1] = dip_1
        ################################################################################

        profl = {}
        profl['e_density']  = field_profile
        profl['glats']      = field_lats
        profl['glons']      = field_lons
        profl['alts']       = field_alts
        profl['ranges']     = np.array(ranges)
        profl['edensTHT']   = edensTHT
        profl['dip']        = profl_dip

        self.profiles[dict_key] = profl 
        return profl

    def generate_wave(self,wave_list=None,sTime=None,currTime=None):
        """
        sTime: Set to self.iri_date if None.
        """

        if wave_list is None:
            wave_list = []
            wave_list.append(dict(src_lat=60.,src_lon=-100.,modulation=0.50,lambda_h=250,T_minutes=15))
            wave_list.append(dict(src_lat=60.,src_lon=-80.,modulation=0.50,lambda_h=250,T_minutes=15))

        rng_from_src = np.zeros([self.lats.size,self.lons.size,self.alts.size])

        lats = np.array(self.lats)
        lats.shape = (lats.size,1,1)

        lons = np.array(self.lons)
        lons.shape = (1,lons.size,1)

        alts = np.array(self.alts)
        alts.shape = (1,1,alts.size)

        edens_0 = self.edens.copy()
        for wave in wave_list:
            src_lat     = wave['src_lat']
            src_lon     = wave['src_lon']
            modulation  = wave['modulation']
            lambda_h    = wave['lambda_h']
            T_minutes   = wave['T_minutes']

            src_lats = np.ones_like(lats) * src_lat
            src_lons = np.ones_like(lons) * src_lon

            rng_from_src = (utils.Re + alts) * utils.geoPack.greatCircleDist(src_lats,src_lons,lats,lons)

            if sTime is None: sTime = self.iri_date

            if (T_minutes is None) or (currTime is None):
                omega = 0.
                t = 0.
            else:
                omega = (2.*np.pi)/(T_minutes*60.)
                t = (currTime - sTime).total_seconds()
                self.date = currTime

            k_h         = 2.*np.pi/lambda_h
            wave        = edens_0 * modulation*np.sin(k_h*rng_from_src - omega*t)
            self.edens  = self.edens + wave

        self.wave_list = wave_list

    def save_profiles(self,keys=None,output_dir='data',filename=None):
        """
        Save density profiles to binary files for use by the raydarn 
        raytracing routine.
        
        Keyword arguments:
        * keys: List of keys of the profiles to be saved.  If None, all profiles
            attached to this object will beb saved.
        * filename: If none, the defaulte fname_base saved in the profile dictionary
            will be used.
        """
        if keys is None: keys = self.profiles.keys()

        paths = []
        for key in keys:
            profl = self.profiles[key]
            date_code = profl['date'].strftime('%Y%m%d_%H%MUT')
            fname_base = '{date_code}_{radar}_{beam:02d}'.format(
                    date_code=date_code,radar=profl['radar'],beam=profl['beam'])

#            Raydarn wants to see the following variables:
#                real*4,dimension(500,500),intent(out)::     edensARR
#                real*4,dimension(500,2),intent(out)::       edensPOS
#                real*4,dimension(500,2),intent(out)::       dip
#                real*4,dimension(500),intent(out)::         edensTHT

            out_array   = np.zeros([500,505],dtype=np.float32)
            out_array[:,:500]       = profl['e_density']
            out_array[:,500]        = profl['glats']
            out_array[:,501]        = profl['glons']
            out_array[:,502:504]    = profl['dip']
            out_array[:,504]        = profl['edensTHT']
            
            if filename is None:
                _filename = os.path.join(output_dir,'{0}.dat'.format(profl['fname_base']))
            else:
                _filename = os.path.join(output_dir,filename)

            # Output binary array in row-major file.
            out_array.T.tofile(_filename)
            paths.append(_filename)

        return paths

    def save_altitude(self,alt=250,output_dir='data',filename=None):
        """
        Save an altitude cut to a pickle file for later plotting.
        
        Keyword arguments:
        * filename: If none, the defaulte fname_base saved in the profile dictionary
            will be used.
        """

        alt_inx     = np.argmin(np.abs(self.alts-alt))
        real_alt    = self.alts[alt_inx]

        dct = {}
        dct['date']     = self.date
        dct['iri_date'] = self.iri_date
        dct['lats']     = self.lats
        dct['lons']     = self.lons
        dct['alt']      = real_alt
        dct['edens']    = self.edens[:,:,alt_inx]

        if hasattr(self,'wave_list'):
            dct['wave_list'] = self.wave_list

        if filename is None:
            fname = '{0}_{1:03.0f}km_edens.p'.format(self.date.strftime('%Y%m%d_%H%MUT'),float(real_alt))
            _filename = os.path.join(output_dir,fname)
        else:
            _filename = os.path.join(output_dir,filename)

        with open(_filename,'wb') as fl:
            pickle.dump(dct,fl)

        return _filename

    def generate_radar_profiles(self,radar,fov=None,site=None,date=None,beams=None,ranges_tup=(0,2500,5)):
        """
        Generate electron density profiles for a radar field of view (FOV).

        radar:      3-letter radar code
        date:       If None, uses self.date.
        beams:      If None, uses all radar beams.
        ranges_tup: (min_range, max_range, range_step) [km]
        """

        ranges = np.arange(*ranges_tup)

        if date is None: date = self.date

        if (fov is None) or (site is None):
            print('WARNING: fov or site object set to None.  Using default field-of-view for {0}.'.format(radar.upper()))
            site      = pydarn.radar.radStruct.site(code=radar,dt=date)
        #   fov       = pydarn.radar.radFov.fov(frang=myBeam.prm.frang, rsep=myBeam.prm.rsep, site=site,elevation=fovElevation,model=fovModel,coords=fovCoords)
            fov       = pydarn.radar.radFov.fov(site=site,model='IS',coords='geo')

        if beams is None: beams = fov.beams # Use all beams from radar.

        for beam in beams:
            az = site.beamToAzim(beam)
            flat_flon = []
            for x in ranges:
                tmp = utils.geoPack.calcDistPnt(site.geolat, site.geolon,origAlt=0,
                        az=az,dist=x,el=0)
                flat_flon.append([tmp['distLat'],tmp['distLon']])

            flat_flon = np.array(flat_flon)

            shape = flat_flon.shape
            flat_flon = flat_flon.flatten()
            tf = flat_flon < 0
            flat_flon[tf] = 360. + flat_flon[tf]
            flat_flon.shape = shape

            dict_key = '{0}-{1:d}'.format(radar,beam)

            profl = self.generate_profile(flat_flon[:,0],flat_flon[:,1],ranges,dict_key)

            # Base filename to be used with this profile.
            date_code   = self.date.strftime('%Y%m%d_%H%MUT')
            fname_base  = '{date_code}_{radar}_{beam:02d}'.format(
                    date_code=date_code,radar=radar,beam=beam)

            # Add extra meta-data.
            profl['date']   = self.date
            profl['radar']  = radar
            profl['beam']   = beam
            profl['fname_base'] = fname_base

        return self.profiles

def rt_fan_plot(rt_run_data,fov=None,figsize=(10,8),
        output_dir='data',filename=None,
        basemap_dict =  {'lat_0':60., 'lon_0':-105, 'width':6500e3, 'height':6500e3, 'lat_ts':50., 'projection':'stere','resolution':'l'},
        multiproc=False,
        nprocs=None):

    fit_df = rt_run_data['fit_df']
    
    ################################################################################
    # Plot fan plot.
    radars  = fit_df.radar.unique()
    dcts    = []
    for radar in radars:
        fit_df_1 = fit_df[fit_df.radar == radar]
        rt_times = fit_df_1.rt_time.unique()
        for rt_time in rt_times:
            rt_time = pd.to_datetime(rt_time).to_datetime()
            fit_df_2 = fit_df_1[fit_df_1.rt_time == rt_time]

            dct = {}
            dct['rt_time']      = rt_time
            dct['fit_df_2']     = fit_df_2
            dct['radar']        = radar
            dct['fov']          = fov
            dct['figsize']      = figsize
            dct['filename']     = filename
            dct['output_dir']   = output_dir
            dct['basemap_dict'] = basemap_dict
            dcts.append(dct)

    if multiproc:
        pool    = multiprocessing.Pool(int(nprocs))
        pool.map(rt_fan_plot_dct,dcts)
        pool.close()
        pool.join()
    else:
        for dct in dcts:
            rt_fan_plot_dct(dct)

def rt_fan_plot_dct(dct):
    return rt_fan_plot_single(**dct)

def rt_fan_plot_single(
        rt_time,
        fit_df_2,
        radar,
        fov,
        figsize,
        filename,
        output_dir,
        basemap_dict):

    fig  = plt.figure(figsize=figsize)
    ax   = fig.add_subplot(111)
    m = Basemap(ax=ax,**basemap_dict)

    # draw parallels and meridians.
#                                                              [Left,  Right, Top,   Bottom]
    m.drawparallels(np.arange(-80.,81.,20.),  color='k',labels=[True,  False, False, False])
    m.drawmeridians(np.arange(-180.,181.,20.),color='k',labels=[False, False, False, True])
    m.drawcoastlines(color='1.',linewidth=4)
    m.drawcoastlines(linewidth=1)
    m.drawmapboundary(fill_color='w')

    if fov is None:
        site    = pydarn.radar.radStruct.site(code=radar,dt=rt_time)
        fov     = pydarn.radar.radFov.fov(site=site,model='IS',coords='geo')

    data_arr    = np.zeros([fov.beams.size,fov.gates.size],dtype=np.float)
    for inx,row in fit_df_2.iterrows():
        beam = row['bmnum']
        beam_inx = np.where(fov.beams == beam)[0][0]
        for gate_inx,pwr in enumerate(row.rt_power_linear):
            if beam_inx == len(fov.beams) or \
               gate_inx == len(fov.gates): continue
            data_arr[beam_inx,gate_inx] = pwr

    log = True
    if log:
        top_val = 30 

        data_arr = 10*np.log10(data_arr/data_arr.max()) + top_val

        shape           = data_arr.shape
        data_arr        = data_arr.flatten()
        tf              = data_arr < 0
        data_arr[tf]    = np.nan
        data_arr.shape  = shape
        
        scale,ticks,cmap = calculate_scale(data_arr,lim=(0,top_val))
        tick_step = 10
        ticks   = np.arange(0,top_val+tick_step,tick_step)
        bounds  = np.linspace(scale[0],scale[1],256)
        norm    = matplotlib.colors.BoundaryNorm(bounds,cmap.N)
    else:
        scale,ticks,cmap = calculate_scale(data_arr)
        bounds  = np.linspace(scale[0],scale[1],256)
        norm    = matplotlib.colors.BoundaryNorm(bounds,cmap.N)

    new_plot = False
    if not new_plot:
        #######################
        #Overlay radar field of view.
        verts = []
        data  = []
        for beam_inx in range(data_arr.shape[0]):
            for gate_inx in range(data_arr.shape[1]):
                if beam_inx == len(fov.beams) or \
                   gate_inx == len(fov.gates): continue

                if not np.isfinite(data_arr[beam_inx,gate_inx]): continue

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
                data.append(data_arr[beam_inx,gate_inx])

        data = np.array(data)
        pcoll = PolyCollection(np.array(verts),edgecolors='face',linewidths=0,closed=False,cmap=cmap,norm=norm)
        pcoll.set_array(data)
        ax.add_collection(pcoll,autolim=False)

#                pcoll_1 = PolyCollection(np.array(verts),facecolor=(1,1,1,0),edgecolors='black',linewidths=None,closed=False,cmap=cmap,norm=norm,zorder=200)
##                    pcoll_1.set_array(data)
#                ax.add_collection(pcoll_1,autolim=False)
        
        for gate_inx,gate in enumerate(fov.gates):
            for beam_inx in [fov.beams.min(),fov.beams.max()+1]:
                lats = np.zeros(2)
                lons = np.zeros(2)
                lats[0] = fov.latFull[beam_inx,gate_inx+0]
                lats[1] = fov.latFull[beam_inx,gate_inx+1]
                lons[0] = fov.lonFull[beam_inx,gate_inx+0]
                lons[1] = fov.lonFull[beam_inx,gate_inx+1]

                m.plot(lons,lats,color='black',latlon=True)

        for beam_inx,beam in enumerate(fov.beams):
            for gate_inx in [fov.gates.min(),fov.gates.max()+1]:
                lats = np.zeros(2)
                lons = np.zeros(2)
                lats[0] = fov.latFull[beam_inx+0,gate_inx]
                lats[1] = fov.latFull[beam_inx+1,gate_inx]
                lons[0] = fov.lonFull[beam_inx+0,gate_inx]
                lons[1] = fov.lonFull[beam_inx+1,gate_inx]

                m.plot(lons,lats,color='black',latlon=True)

    else:
        LONS    = fov.lonFull
        LATS    = fov.latFull
        pcoll = m.pcolormesh(LONS,LATS,data_arr,cmap=cmap,norm=norm,zorder=200,latlon=True,edgecolors='black')

    cbar    = fig.colorbar(pcoll,orientation='vertical',shrink=0.60,pad= 0.015,ticks=ticks)
    txt     = 'Ray Trace GS Power [dB]'
    cbar.set_label(txt)

    m.nightshade(rt_time)
    ax.set_title('{0}'.format(radar.upper()))

    freqs = fit_df_2['prm_tfreq'].unique()/1000.
    if len(freqs) == 1:
        freq_str = '{:.03f} MHz'.format(freqs[0])
    else:
        freq_str = '{!s} MHz'.format(freqs)

    txt = []
    txt.append('SuperDARN Ray Trace Simulation')
    txt.append('{} ({}) {}'.format(radar.upper(),freq_str,rt_time.strftime('%d %b %Y %H%M UT')))
    ax.set_title('\n'.join(txt))

    fig.tight_layout()

    if filename is None:
        fname = '{0}_{1}_rt_map.png'.format(rt_time.strftime('%Y%m%d_%H%MUT'),radar)
        _filename = os.path.join(output_dir,fname)
    else:
        _filename = os.path.join(output_dir,filename)

    fig.savefig(_filename,bbox_inches='tight')
    plt.close()

def df_log_power(rt_run,fov,noise_floor=30.,set_nans=True):
    """
    Compute the log power in an rt_run dataframe.
    
    rt_run:     dictionary object created by the MSTID raytracer script.
    fov:        Field-of-view object for the radar.  Assume only one radar in
                the rt_run object for now.
    noise_floor: dB level to set as highest value.
    set_nans:   If True, set everything below 0 dB as np.nan.
    """

    fit_df  = rt_run['fit_df']

    # Find the dates and make sure they are actually datetime.datetime objects.
    rt_times   = fit_df['rt_time'].unique()
    rt_times   = [pd.to_datetime(x).to_datetime() for x in rt_times]
    rt_times.sort()
    rt_times   = np.array(rt_times)

    nbeams  = fov.beams.size
    ngates  = fov.gates.size

    # Create numpy array to hold power data.
    pwr_arr = np.zeros([nbeams,ngates,rt_times.size],dtype=np.float)

    # Don't assume things are in order... explicitly identify which index each thing
    # belongs to.
    df_inxs = np.zeros([nbeams,rt_times.size],dtype=np.int)    # Keep track of what row in the data
                                                            # frame each iteration came from.
    for df_inx,row in fit_df.iterrows():
        bmnum_inx = row['bmnum']
        rt_time_inx = np.where(row['rt_time'] == rt_times)[0][0]

        df_inxs[bmnum_inx,rt_time_inx]   = df_inx
        pwr_arr[bmnum_inx,:,rt_time_inx] = row['rt_power_linear']

    ########## 
    # Now manipulate the data.
    log_pwr = 10*np.log10(pwr_arr/pwr_arr.max()) + noise_floor
    if set_nans:
        shape           = log_pwr.shape
        log_pwr         = log_pwr.flatten()
        tf              = log_pwr < 0
        log_pwr[tf]     = np.nan
        log_pwr.shape   = shape

    ########## 
    # Put everything back into the data frame and return it.
    fit_df['rt_power_log']   = [np.zeros(ngates,dtype=np.float)]*len(fit_df)

    for inx in fit_df.index:
        bmnum_inx,rt_time_inx = np.where(inx == df_inxs)
        fit_df.set_value(inx,'rt_power_log',log_pwr[bmnum_inx[0],:,rt_time_inx[0]])

    return rt_run

class BeamList(object):
    def __init__(self,rt_run_data,sTime=None,eTime=None):
        fit_df = rt_run_data['fit_df']

        if sTime is not None:
            fit_df = fit_df[fit_df.rt_time >= sTime]

        if eTime is not None:
            fit_df = fit_df[fit_df.rt_time < eTime]

        rt_times = fit_df.rt_time.unique()
        rt_times = [pd.to_datetime(x).to_datetime() for x in rt_times]
        rt_times.sort()
        rt_times = np.array(rt_times)

        self.sTime = rt_times.min()
        self.eTime = rt_times.max()

        self.radar = fit_df.iloc[0]['radar']
        self.stid  = fit_df.iloc[0]['prm_stid']
        self.cp    = fit_df.iloc[0]['prm_cpid']
        self.channel = 1
        self.dType = 'raytrace'
        self.fType = 'raytrace'


        beam_list = [] 
        for inx,row in fit_df.iterrows():
            myBeam = pydarn.sdio.beamData()

            myBeam.bmnum        = row['bmnum']   # Beam number
            myBeam.fPtr         = None          # File pointer (not used)
            myBeam.fitex        = None          # Fitex flag
        #    myBeam.fit          = object 
        #    myBeam.prm          = object 
            myBeam.recordDict   = None
            myBeam.stid         = row['prm_stid']
            myBeam.lmfit        = None 
            myBeam.exflg        = None 
            myBeam.iqflg        = None 
            myBeam.offset       = None
        #    myBeam.rawacf       = object 
            myBeam.lmflg        = None 
            myBeam.rawflg       = None 
            myBeam.fType        = 'fitacf'
            myBeam.time         = row['rt_time']
            myBeam.acflg        = None 
            myBeam.cp           = row['prm_cpid']
        #    myBeam.iqdat        = object 
            myBeam.fitacf       = None 
            myBeam.channel      = self.channel

            # PRM Block #################################################################### 
            myBeam.prm.bmazm        = row['prm_bmazm']
            myBeam.prm.lagfr        = row['prm_lagfr']
            myBeam.prm.nave         = row['prm_nave']
            myBeam.prm.rsep         = row['prm_rsep']
            myBeam.prm.noisemean    = row['prm_noise']['mean']
            myBeam.prm.rxrise       = row['prm_rxrise']
            myBeam.prm.xcf          = row['prm_xcf']
            myBeam.prm.frang        = row['prm_frang']
            myBeam.prm.mpinc        = row['prm_mpinc']
            myBeam.prm.noisesearch  = row['prm_noise']['search']
            myBeam.prm.scan         = np.where(row['rt_time'] == rt_times)[0][0]
            myBeam.prm.ifmode       = row['prm_ifmode']
            myBeam.prm.mplgexs      = row['prm_mplgexs']
            myBeam.prm.smsep        = row['prm_smsep']
            myBeam.prm.inttsc       = row['prm_intt']['sc']
            myBeam.prm.mplgs        = row['prm_mplgs']
            myBeam.prm.nrang        = row['prm_nrang']
            myBeam.prm.tfreq        = row['prm_tfreq']
            myBeam.prm.inttus       = row['prm_intt']['us']
            myBeam.prm.mppul        = row['prm_mppul']
            myBeam.prm.txpl         = row['prm_txpl']
            myBeam.prm.ptab         = row['prm_pulse']
            myBeam.prm.ltab         = row['prm_lag']

            myBeam.prm.noisesky     = row['fit_noise']['sky']

            # FIT Block ####################################################################

            tf      = np.isfinite(row['rt_power_log'])
            slist   = np.where(tf)[0]
            p_l_arr = np.array(row['rt_power_log'])[slist]
            zeros   = [0.]*p_l_arr.size

            myBeam.fit.p_l          = p_l_arr.tolist()
            myBeam.fit.phi0_e       = zeros 
            myBeam.fit.v            = zeros
            myBeam.fit.w_s_e        = zeros
            myBeam.fit.elv          = np.array(row['rt_elv'])[slist].tolist()
            myBeam.fit.p_l_e        = zeros
            myBeam.fit.pwr0         = np.array(row['rt_pwr0'])[slist].tolist()
            myBeam.fit.v_e          = zeros
            myBeam.fit.gflg         = np.array(row['rt_gflg'])[slist].tolist()
            myBeam.fit.p_s          = zeros
            myBeam.fit.qflg         = np.array(row['rt_qflg'])[slist].tolist()
            myBeam.fit.w_l          = zeros
            myBeam.fit.nlag         = row['fit_nlag']
            myBeam.fit.p_s_e        = zeros
            myBeam.fit.slist        = slist.tolist()
            myBeam.fit.w_l_e        = zeros
            myBeam.fit.npnts        = None
            myBeam.fit.phi0         = zeros
            myBeam.fit.w_s          = zeros
            
            beam_list.append(myBeam)

        self.beam_list  = beam_list
        self.scan_index = 0

    def readScan(self):
        scan = []
        for item in self.beam_list:
            if item.prm.scan == self.scan_index:
                scan.append(item)

        if scan == []:
            scan = None
        else:
            self.scan_index += 1

        return scan

    def close(self):
        """ This method is here for compatibility with myPtr(). """
        pass

def periodic_extend(rt_run_data,extend_sTime,extend_eTime):
    fit_df = rt_run_data['fit_df']
    
    new_dt  = (extend_eTime - extend_sTime).total_seconds()
    old_dt  = (fit_df.rt_time.max() - fit_df.rt_time.min()).total_seconds()
    cycles  = int(new_dt/old_dt)

    rt_times = fit_df.rt_time.unique()
    rt_times = [pd.to_datetime(x).to_datetime() for x in rt_times]
    rt_times.sort()
    rt_times = np.array(rt_times)

    delt    = rt_times[1] - rt_times[0]

    new_list    = []
    time_counter = 0
    for cycle in range(cycles):
        for rt_time_inx,rt_time in enumerate(rt_times):
            curr_time   = extend_sTime + delt*time_counter
            if curr_time >= extend_eTime: break

            t_dict = {'hr': curr_time.hour, 'mt': curr_time.minute, 'mo': curr_time.month, 'yr': curr_time.year, 'dy': curr_time.day}
            scan = fit_df[fit_df.rt_time == rt_time]
            for row_inx,row in scan.iterrows():
                row['rt_time']  = curr_time
                row['prm_time'] = t_dict
                new_list.append(row)
            time_counter += 1

    fit_df = pd.DataFrame(new_list,index=np.arange(len(new_list)))
    rt_run_data['fit_df'] = fit_df
    return rt_run_data
