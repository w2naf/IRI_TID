#!/usr/bin/env python
# coding: utf-8
import sys
import os
import shutil
import datetime
import pickle
import glob
import copy

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt

from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

import subprocess
import multiprocessing

from davitpy import pydarn
from davitpy.models import raydarn
from davitpy import utils
import iono_lib
from iono_lib import iono_3d

# Setup logging. ################################################################
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

fh = logging.FileHandler('log_filename.txt',mode='w')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

def run_mstid_rt(
    radar                                       = 'bks',
    ngates                                      = 55,
    freq                                        = 10.5,
    sTime                                       = datetime.datetime(2014,11,1, 16, 0),
    eTime                                       = datetime.datetime(2014,11,1, 16, 4),
    delta_t                                     = datetime.timedelta(minutes=2),
    extend                                      = False,
    extend_sTime                                = datetime.datetime(2014,11,1),
    extend_eTime                                = datetime.datetime(2014,11,2),
    wave_list                                   = None,
    event_note                                  = None,
    clear_output_dirs                           = True,
    run_get_ionosphere                          = True,
    load_background_ionosphere_from_pickle      = True,
    background_ionosphere_source_for_symlink    = None,
    load_perturbed_ionosphere_from_pickle       = True,
    save_background_ionosphere_to_pickle        = True,
    save_perturbed_ionosphere_to_pickle         = False,
    run_raytrace                                = True,
    plot_rtmaps                                 = True,
    base_dir                                    = 'data',
    iono_dict                                   = None,
    subproc                                     = False,
    multiproc                                   = True,
    nprocs                                      = 6,
    **kwargs
        ):
    """
    Compute a 4D IRI Run and ray trace through it.  Generate a bunch of plots along
    the way, and output data files/objects that can be fed into the MUSIC code
    for further analysis.

    *extend: Repeat the ionosphere computed between sTime and eTime periodically to
             fill from extend_sTime to exentd_eTime.  This is a quick way to generate
             a long duration wave without having to do a long-duration IRI run.
             Of course, you lose the diurnal variations included in a normal IRI run
             when you do this.

    * wave_list: A list of dictionaries defining waves with which to perturb the IRI.
                 The wave is modeled as a circular sin function traveling out from the
                 src_lat,src_lon: wave = edens_0 * modulation*np.sin(k_h*rng_from_src - omega*t)

      The list should be of the form:
      [{'src_lat': 59.9,  'src_lon':  127.8, 'modulation': 0.3, 'T_minutes': 40, 'lambda_h': 451}]

      where:
        src_lat:    Latitude of the wave source.
        src_lon:    Longitude of the wave source.
        modulation: Amplitude of wave.
        T_minutes:  Period in minutes of wave [minutes].
        lambda_h:   Horizontal wavelength [km]

    *event_note: string appended to automatically generated directory names/filenames

    *clear_output_dirs: (bool) Remove any old data before starting

    *load_background_ionosphere_from_pickle: (bool) Use a pre-computed ionosphere from a pickle file.  This
        precomputed ionosphere is normally stored in a non-perturbed state.

    """

    info_dict   = {}    # Dictionary to save things that we might want the
                        # funtion to return.

    ################################################################################
    # Set run parameters. ##########################################################
    ################################################################################
    basemap_dict        =  {'lat_ts':50., 'projection':'stere', 'resolution':'l'}

    # Wide area North America
    basemap_dict['lat_0']   =   60.
#    basemap_dict['lon_0']   = -105.
    basemap_dict['lon_0']   = -70.
    basemap_dict['width']   =  6500e3
    basemap_dict['height']  =  6500e3

#    # Tight Crop US
#    basemap_dict['lat_0']   =  50.
#    basemap_dict['lon_0']   = -90.
#    basemap_dict['width']   =  6500e3
#    basemap_dict['height']  =  4000e3

#    wave_list = []
#    wave_list.append(dict(src_lat=61.64,src_lon=-136.77,modulation=0.30,lambda_h=250,T_minutes=20))  #Ortho at 4000 km
    #wave_list.append(dict(src_lat=55.78,src_lon=-31.12,modulation=0.30,lambda_h=250,T_minutes=20))  #Parallel at 4000 km
    #wave_list.append(dict(src_lat=53.91,src_lon=-54.41,modulation=0.30,lambda_h=250,T_minutes=20))  #Parallel at 2500 km
    #wave_list.append(dict(src_lat=57.42,src_lon=-111.32,modulation=0.30,lambda_h=250,T_minutes=20)) #Ortho at 2500 km
    #wave_list.append(dict(src_lat=60.,src_lon=-100.,modulation=0.50,lambda_h=250,T_minutes=15))
    #wave_list.append(dict(src_lat=45.,src_lon=-80.,modulation=0.50,lambda_h=250,T_minutes=15))

    ################################################################################
    # Initial Set-Up ###############################################################
    ################################################################################

    #event_code  = '{0}_{1}-{2}'.format(radar,sTime.strftime('%Y%m%d_%H%M'),eTime.strftime('%Y%m%d_%H%M'))
    event_code  = '{0}_{1}'.format(radar,sTime.strftime('%Y%m%d'))

    if event_note is not None:
        event_code = '_'.join([event_code,event_note])

    final_dir   = os.path.join(base_dir,event_code)
    dirs        = iono_lib.make_dirs(final_dir,clear_output_dirs=clear_output_dirs)

    # Link to a different background ionosphere run.
    if load_background_ionosphere_from_pickle and background_ionosphere_source_for_symlink is not None:
        dst = dirs['background_edens_raw']

        try:
            if os.path.isdir(dst) and not os.path.islink(dst):
                shutil.rmtree(dst)
            elif os.path.islink(dst):
                os.remove(dst)
        except:
            pass

        src = os.path.abspath(os.path.join(background_ionosphere_source_for_symlink,'background_edens_raw'))
        os.symlink(src,dst)

    site    = pydarn.radar.radStruct.site(code=radar,dt=sTime)
    #fov    = pydarn.radar.radFov.fov(frang=myBeam.prm.frang, rsep=myBeam.prm.rsep, site=site,elevation=fovElevation,model=fovModel,coords=fovCoords)
    frang   = 180.
    rsep    = 45.
    fov     = pydarn.radar.radFov.fov(frang=frang,rsep=rsep,site=site,model='IS',coords='geo',ngates=ngates)

    # Make sure ngates is int
    ngates  = int(ngates)

    time_steps = [sTime]
    while time_steps[-1] < eTime:
        new_time = time_steps[-1] + delta_t
        if new_time >= eTime:
            break
        else:
            time_steps.append(time_steps[-1]+delta_t)
    time_steps = np.array(time_steps)

    # Generate Ionosphere ##########################################################
    if run_get_ionosphere:
        iono_runfiles = []
        for time_step_inx,time_step in enumerate(time_steps):
            plot_perturbed_profiles = True if time_step_inx == 0 else False

            generate_ionosphere_dict    = dict(
                radar           = radar,
                fov             = fov,
                site            = site,
                sTime           = sTime,
                time_step       = time_step,
                dirs            = dirs,
                load_background_ionosphere_from_pickle      = load_background_ionosphere_from_pickle,
                save_background_ionosphere_to_pickle        = save_background_ionosphere_to_pickle,
                iono_dict                                   = iono_dict,
                load_perturbed_ionosphere_from_pickle       = load_perturbed_ionosphere_from_pickle, 
                save_perturbed_ionosphere_to_pickle         = save_perturbed_ionosphere_to_pickle,
                wave_list                                   = wave_list,
                plot_perturbed_profiles                     = plot_perturbed_profiles,
                plot_perturbed_maps                         = True,
                basemap_dict                                = basemap_dict)

            time_str    = time_step.strftime('%Y%m%d_%H%MUT')
            pkl_fname   = '_'.join([time_str,'background_edens_runfile.p'])
            pkl_fpath   = os.path.join(dirs['background_edens_runfiles'],pkl_fname)
            with open(pkl_fpath,'wb') as fl:
                pickle.dump(generate_ionosphere_dict,fl)

            iono_runfiles.append(pkl_fpath)

        if multiproc:
            pool    = multiprocessing.Pool(int(nprocs))
            pool.map(run_generate_ionosphere_runfile,iono_runfiles)
            pool.close()
            pool.join()
        else:
            for pkl_fpath in iono_runfiles:
                run_generate_ionosphere_runfile(pkl_fpath,subproc=subproc)

        if not run_raytrace: return

    # Ray Tracing Section ##########################################################
    if run_raytrace:
        edens_files = glob.glob(os.path.join(dirs['edens_profl_dats'],'*.dat'))
        edens_files.sort()
        
        # Save the input parameters for each Ray Trace run into individual pickle files.
        rtrunfile_paths = []
        fitdict_fpaths   = []
        plot_rt_time    = None
        for edens_inx,edens_file in enumerate(edens_files):
            edens_path      = os.path.abspath(edens_file)
            basename        = os.path.basename(edens_file)
            base_path       = os.path.split(os.path.split(edens_path)[0])[0]
            key             = os.path.splitext(basename)[0]
            rt_data_path    = os.path.join(base_path,'rt_data',key)
            rtrunfile_path  = os.path.join(rt_data_path,'{}_rtrunfile.p'.format(key))
            fitdict_fpath    = os.path.join(rt_data_path,'{}_fitdict.p'.format(key))

            try:
                os.makedirs(rt_data_path)
            except:
                pass

            tmp = {}
            tmp['key']                          = key
            tmp['edens_path']                   = edens_path
            tmp['edens_file']                   = edens_file
            tmp['rt_time']                      = datetime.datetime.strptime(basename[:13],'%Y%m%d_%H%M')
            tmp['radar']                        = basename[16:19]
            tmp['freq']                         = freq
            tmp['beam']                         = int(basename[20:22])
            tmp['ngates']                       = ngates
            tmp['site']                         = site
            tmp['fov']                          = fov
            tmp['frang']                        = frang
            tmp['rsep']                         = rsep
            tmp['rt_data_path']                 = rt_data_path
            tmp['rt_profl_figs']                = os.path.join(base_path,'rt_profl_figs')

            # Only plot profiles for the first time step of the run.
            if plot_rt_time is None:
                plot_rt_time = tmp['rt_time']

            if tmp['rt_time'] == plot_rt_time:
                tmp['plot_profiles'] = True

            with open(rtrunfile_path,'wb') as fl:
                pickle.dump(tmp,fl)

            rtrunfile_paths.append(rtrunfile_path)
            fitdict_fpaths.append(fitdict_fpath)

        # Run the ray trace and save output data to disk.
        if multiproc and subproc:
            pool    = multiprocessing.Pool(int(nprocs))
            pool.map(run_ray_trace_edens_profl_runfile,rtrunfile_paths)
            pool.close()
            pool.join()
        else:
            for rtrunfile_path in rtrunfile_paths:
                run_ray_trace_edens_profl_runfile(rtrunfile_path,subproc=subproc)

        # Aggregate all raytrace results into a single file that can has the equivalent contents of a fit
        # file (at least in terms of ground power scatter).

        fit_list    = []
        for fitdict_fpath in fitdict_fpaths:
            with open(fitdict_fpath,'rb') as fl:
                fitdict = pickle.load(fl)
            fit_list.append(fitdict)

        fit_df = pd.DataFrame(fit_list)
        rt_run_data = {}
        rt_run_data['fit_df']       = fit_df
        rt_run_data['wave_list']    = wave_list

        filename    = os.path.join(dirs['000_general_output'],event_code+'_rtresults.p')
        with open(filename,'wb') as fl:
            pickle.dump(rt_run_data,fl)
    else:
        filename    = os.path.join(dirs['000_general_output'],event_code+'_rtresults.p')
        with open(filename,'rb') as fl:
            rt_run_data = pickle.load(fl)

    if plot_rtmaps:
        logger.debug('Plotting ray trace fan plots.')
        iono_lib.rt_fan_plot(rt_run_data,fov=fov,basemap_dict=basemap_dict,output_dir=dirs['rt_maps'],
                multiproc=multiproc,nprocs=nprocs)

    if extend:
        logger.debug('Creating a periodic expansion from {0} to {1}.'.format(str(extend_sTime),str(extend_eTime)))
        rt_run_data = iono_lib.periodic_extend(rt_run_data,extend_sTime,extend_eTime)

    logger.debug('Converting linear power to log power.')
    rt_run_data = iono_lib.df_log_power(rt_run_data,fov)

    filename    = os.path.join(dirs['000_general_output'],event_code+'_logp_rtresults.p')
    with open(filename,'wb') as fl:
        pickle.dump(rt_run_data,fl)

    logger.debug('Forming dataframe into beam_list and saving to pickle.')
    beam_list   = iono_lib.BeamList(rt_run_data)
    filename    = os.path.join(dirs['000_general_output'],'{0}_beam_list.p'.format(event_code))
    info_dict['beam_list_fname']    = filename
    with open(filename,'wb') as fl:
        pickle.dump(beam_list,fl)

    logger.debug('Converting to MUSIC object and saving to pickle.')
    dataObj     = pydarn.proc.music.musicArray(beam_list,fovModel='GS')
    filename    = os.path.join(dirs['000_general_output'],'{0}_music_obj.p'.format(event_code))
    with open(filename,'wb') as fl:
        pickle.dump(dataObj,fl)

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    pydarn.plotting.musicFan(dataObj,axis=ax)
    filename    = os.path.join(dirs['000_general_output'],'map.png')
    fig.savefig(filename)
    plt.close()
    
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    pydarn.plotting.musicRTI(dataObj,axis=ax)
    filename    = os.path.join(dirs['000_general_output'],'rti.png')
    fig.savefig(filename)
    plt.close()

    return info_dict

#    shutil.move(dirs['base'],os.path.join('data',event_code))

def get_source_loc(rng,azm,radar,sTime,lat_ctr=None,lon_ctr=None,azm_plus_boresite=False):
    """
    azm is defined as where the MSTID is going; not where it is coming from.
    """
    site    = pydarn.radar.radStruct.site(code=radar,dt=sTime)
    fov     = pydarn.radar.radFov.fov(site=site,model='IS',coords='geo')

    if lat_ctr is None and lon_ctr is None:
        if radar == 'bks':
            lat_ctr =  42.8
            lon_ctr = -84.2

        else:
            print 'You must specify a radar center lat/lon!'
            import ipdb; ipdb.set_trace()

    if azm_plus_boresite:
        azm = azm + site.boresite

    lat,lon = utils.greatCircleMove(lat_ctr,lon_ctr,rng,azm+180.)
    return lat,lon,azm

def generate_ionosphere(radar,fov,site,sTime,time_step,dirs,
        load_background_ionosphere_from_pickle  	= False,
        save_background_ionosphere_to_pickle    	= True,
        iono_dict                               	= None,
        load_perturbed_ionosphere_from_pickle   	= False,
        save_perturbed_ionosphere_to_pickle     	= False,
        wave_list                               	= None,
        plot_perturbed_profiles                 	= True,
        plot_perturbed_maps                     	= True,
        basemap_dict                            	= None):

    ##### Load background ionosphere.
    time_str    = time_step.strftime('%Y%m%d_%H%MUT')
    pkl_fname   = '_'.join([time_str,'background_edens_raw.p'])
    pkl_fpath   = os.path.join(dirs['background_edens_raw'],pkl_fname)

    t0 = datetime.datetime.now()
    logger.debug('Getting RAW IONOSPHERE for {}...'.format(time_str))
    if load_background_ionosphere_from_pickle and os.path.exists(pkl_fpath):
        logger.debug('Loading: {}'.format(pkl_fpath))
        with open(pkl_fpath,'rb') as fl:
            iono_0 = pickle.load(fl)
    elif not os.path.exists(pkl_fpath):
        if iono_dict is None:
            iono_dict = dict(
                    hgt_0 =   60., hgt_1 =  560., hgt_step=1.0,
                    lat_0 =   32., lat_1 =   80., lat_step=1.00,
                    lon_0 = -100., lon_1 =  -40., lon_step=1.00)

        iono_0          = iono_lib.iono_3d(time_step,**iono_dict)
        iono_0.iri_date = iono_0.date

        if save_background_ionosphere_to_pickle:
            logger.debug('Saving: {}'.format(pkl_fpath))
            with open(pkl_fpath,'wb') as fl:
                pickle.dump(iono_0,fl)
    t1 = datetime.datetime.now()
    logger.debug('Total RAW IONOSPHERE get time: {0} s'.format(str(t1-t0)) )

    ##### Generate waves in ionosphere.
    time_str    = time_step.strftime('%Y%m%d_%H%MUT')
    pkl_fname   = '_'.join([time_str,'background_edens_perturbed.p'])
    pkl_fpath   = os.path.join(dirs['background_edens_perturbed'],pkl_fname)
    t0 = datetime.datetime.now()

    logger.debug('Getting PERTURBED IONOSPHERE for {}...'.format(time_str))
    if load_perturbed_ionosphere_from_pickle and os.path.exists(pkl_fpath):
        del iono_0
        logger.debug('Loading: {}'.format(pkl_fpath))
        with open(pkl_fpath,'rb') as fl:
            iono = pickle.load(fl)
    else:
        logger.debug('Generating ionospheric profiles for {0}'.format(str(time_step)))
        iono = iono_0
        if wave_list is not None:
            iono.generate_wave(wave_list=wave_list,sTime=sTime,currTime=time_step)
        iono.generate_radar_profiles(radar,fov=fov,site=site)
        iono.save_profiles(output_dir=dirs['edens_profl_dats'])
        iono.save_altitude(output_dir=dirs['edens_alt_p'])

        if save_perturbed_ionosphere_to_pickle:
            logger.debug('Saving: {}'.format(pkl_fpath))
            with open(pkl_fpath,'wb') as fl:
                pickle.dump(iono,fl)
    t1 = datetime.datetime.now()
    logger.debug('Total PERTURBED IONOSPHERE get time: {0} s'.format(str(t1-t0)) )

    ##### Plot perturbed ionosphere.
    if plot_perturbed_profiles:
        t0 = datetime.datetime.now()
        logger.debug('Plotting profiles: {!s}'.format(time_step))
        iono.plot_profiles(output_dir=dirs['edens_profl_figs'])
        t1 = datetime.datetime.now()
        logger.debug('Total profile plotting time: {0} s'.format(str(t1-t0)) )

    if plot_perturbed_maps:
        logger.debug('Plotting ionosphere map.')
        iono.plot_map(radar=radar,fov=fov,output_dir=dirs['edens_maps'],basemap_dict=basemap_dict)

def run_generate_ionosphere_runfile(pkl_fpath,subproc=True):
    """
    Script to launch ionospheric generator as its own subprocess, to prevent
    memory leaks associated with working with large pickle files on python 2
    in linux.
    """
    if subproc:
        cmd = ['./run_generate_ionospere_subprocess.py',pkl_fpath]
        print ' '.join(cmd)
        subprocess.check_call(cmd)
    else:
        with open(pkl_fpath,'rb') as fl:
            generate_ionosphere_dict = pickle.load(fl)

        generate_ionosphere(**generate_ionosphere_dict)

def run_ray_trace_edens_profl_runfile(pkl_fpath,subproc=True):
    """
    Script to raytracer as its own subprocess, to prevent
    memory leaks associated with working with large pickle files on python 2
    in linux.
    """
    if subproc:
        cmd = ['./run_ray_trace_edens_profl_subprocess.py',pkl_fpath]
        print ' '.join(cmd)
        subprocess.check_call(cmd)
    else:
        with open(pkl_fpath,'rb') as fl:
            rt_runfile_dict = pickle.load(fl)

        ray_trace_edens_profl(**rt_runfile_dict)


def ray_trace_edens_profl(
    key,
    edens_path,
    edens_file,
    rt_time,
    radar,
    freq,   # in MHz
    beam,
    ngates,
    site,
    fov,
    frang,
    rsep,
    rt_data_path,
    rt_profl_figs,
    plot_profiles               = False):

    rto_fpath       = os.path.join(rt_data_path,'{}_rto.p'.format(key))
    fitdict_fpath    = os.path.join(rt_data_path,'{}_fitdict.p'.format(key))

    # Do not run computation if we already have a fitdict computed.
    if os.path.exists(fitdict_fpath):
        return

    logger.debug('Now Ray Tracing: {0}'.format(key))
    shutil.rmtree(rt_data_path,ignore_errors=True)
    os.makedirs(rt_data_path)

    # Run the ray tracing for the specified period, radar, beam and frequency
    rto = raydarn.RtRun(sTime=rt_time, eTime=rt_time, rCode=radar, beam=beam, 
            freq=freq, outDir=rt_data_path, nprocs=1,edens_file=edens_file)
    with open(rto_fpath,'wb') as fl:
        pickle.dump(rto,fl)

    # Read rays, electron densities, and scatter into memory
    rto.readRays()
    rto.readEdens()
    rto.readScatter()

    power   = rto.scatter.gate_scatter(beam=beam,fov=fov)
    rt_gflg = [int(x) for x in (power > 0)]

    # Create parameter block. #### 
    # Using parameter block values based on Sebastien's IDL routine rt_makefit.pro
    # Some of the numbers may need to be reviewed.
    #-----------------------------------------
    #- Fill structure: PRM
    #-----------------------------------------
    prm = {}
    prm['prm_revision'] = {}
    prm['prm_revision']['major'] = 2
    prm['prm_revision']['minor'] = 0

    prm['prm_origin'] = {}
    prm['prm_origin']['time'] = datetime.datetime.now()
    prm['prm_origin']['command'] = 'rt_makefit_naf'

    prm['prm_stid']     = pydarn.radar.network().getRadarByCode(radar).id

    prm['prm_time'] = {}
    prm['prm_time']['yr'] = rt_time.year
    prm['prm_time']['mo'] = rt_time.month
    prm['prm_time']['dy'] = rt_time.day
    prm['prm_time']['hr'] = rt_time.hour
    prm['prm_time']['mt'] = rt_time.minute

    prm['prm_lagfr'] = ((frang*1e3*2.)/3.e8)*1e6 # Lag time to first range gate in microsec.
    prm['prm_smsep'] = ((rsep *1e3*2.)/3.e8)*1e6 # Separation between gates in microsec.
    prm['prm_noise'] = {}
    prm['prm_noise']['search'] = 1e-6
    prm['prm_noise']['mean'] = 1e-6
    prm['prm_bmazm'] = site.beamToAzim(beam)
    prm['prm_rxrise'] = site.recrise
    prm['prm_ifmode'] = None    # Possible values --> -1?

    prm['prm_intt'] = {}
    prm['prm_intt']['sc'] = 5
    prm['prm_intt']['us'] = 318181

    prm['prm_txpl']  = 300  # Tranmitted powerl level?
    prm['prm_mplgexs'] = 0 # ?
    prm['prm_mpinc'] = 1500 # ?
    prm['prm_mppul'] = 8    # ?
    prm['prm_mplgs'] = 23   # ?

    prm['prm_nave'] = 1

    prm['prm_nrang'] = fov.gates.size
    prm['prm_frang'] = frang
    prm['prm_rsep']  = rsep
    prm['prm_xcf']   = 1
    prm['prm_tfreq'] = freq*1e3
    prm['prm_mxpwr'] = 1073741824   # ?
    prm['prm_lvmax'] = 20000        # ?

#        prm_pulse = [ 0, 14, 22, 24, 27, 31, 42, 43]    # Pulse table
    prm_pulse = [0]*8   # Pulse table
    prm['prm_pulse'] = prm_pulse

    prm_lag = np.zeros([24,2],dtype=np.int) # Lag table
#        prm_lag[:,0] = np.array([0, 42, 22, 24, 27, 22, 24, 14, 22, 14, 31, 31, 14,  0, 27, 27, 14, 24, 24, 22, 22,  0,  0, 43])
#        prm_lag[:,1] = np.array([0, 43, 24, 27, 31, 27, 31, 22, 31, 24, 42, 43, 27, 14, 42, 43, 31, 42, 43, 42, 43, 22, 24, 43])
    prm['prm_lag'] = prm_lag.tolist()
    
    prm['prm_cpid'] = -7373

    # End fill structure PRM
    #-----------------------------------------

    #-----------------------------------------
    #- Fill structure: FIT
    #-----------------------------------------
    fit = {}
    fit['fit_revision'] = {}
    fit['fit_revision']['major'] = 2
    fit['fit_revision']['minor'] = 0

    fit['fit_noise'] = {}
    fit['fit_noise']['sky']     = 1e-6
    fit['fit_noise']['lag0']    = 1e-6
    fit['fit_noise']['vel']     = 1e-6
    fit['fit_nlag'] = [0]*5     # Sebastien's vals: [6, 7, 17, 8, 5]

    # End fill structure FIT
    #-----------------------------------------

    rt_results = {}
    rt_results['radar']             = radar
    rt_results['bmnum']             = beam
    rt_results['rt_time']           = rt_time
    rt_results['rt_power_linear']   = power
    rt_results['rt_gflg']           = rt_gflg
     
    rt_results['rt_pwr0']           = [0]*ngates    #lagpower
    rt_results['rt_qflg']           = rt_gflg

    rt_results['rt_x_gflg']         = rt_gflg
    rt_results['rt_x_qflg']         = rt_gflg

    rt_results['rt_elv']            = [0]*ngates
    rt_results['rt_elv_low']        = [0]*ngates
    rt_results['rt_elv_high']       = [0]*ngates

    # Save fitdict to disk.
    fitdict         = dict(prm.items() + fit.items() + rt_results.items())
    with open(fitdict_fpath,'wb') as fl:
        pickle.dump(fitdict,fl)

    # Only plot ray trace ionospheric profiles for one scan of a given radar.
    if plot_profiles:
        logger.debug('Plotting ray trace profile: {0}'.format(key))
        # Plot rays and electron densities together
        # Plot range markers (every 250 km)
        fig = plt.figure(figsize=(15,5))
        plt.rcParams.update({'font.size': 14})
        ax, aax, cbax = rto.ionos.plot(rt_time)
        ax, aax, cbax = rto.rays.plot(rt_time, step=10, ax=ax, aax=aax)
        ax.grid()
        time_str    = rt_time.strftime('%d %b %Y %H%M UT')
        txt         = 'Raytrace: {} Beam {:.0f} [{}]'.format(radar.upper(),beam,time_str)
        ax.set_title(txt)
        fname = os.path.join(rt_profl_figs,'{0}_rt_profl.png'.format(key))
        fig.savefig(fname)
        plt.close('all')
