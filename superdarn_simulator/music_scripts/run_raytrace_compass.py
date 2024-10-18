#!/usr/bin/env python
# coding: utf-8

# This is a special version the MSTID raytrace simulator that:
#   a. Will use a pre-calculated background ionophere from a different run by creating 
#       simlinks to it.
#   b. Will run a wave for every cardinal and subcardinal direction to test radar sensitivity.

# General instructions for running the MSTID raytrace simulator.
# 1. Run this script with the desired radar/wave parameters.  This script will:
#   a. Run the IRI and populate an electron density array.
#   b. Perturb the IRI run with the appropriate wave.
#   c. Plot profiles and maps of the IRI run.
#   d. Ray-trace the selected radar through the perturbed IRI run.
#   e. Estimate the ground scatter power seen by the radar.
#   f. Plot profiles and maps of the IRI run.
#   g. Create lots of pickle files along the way of all the calculated results.
#      Importantly, it generates a *_beam_list.p file, which essentially contains
#      the ground scatter information of a fitacf file.  You can pass this to the
#      run_music_on_raytrace_data.py routine as a srcPath, which will then allow the
#      simulation data to be run through the MUSIC processing scheme.

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

from davitpy import pydarn
from davitpy.models import raydarn
from davitpy import utils

from mstid_raytrace import run_mstid_rt, get_source_loc

if __name__ == '__main__':
    iono_dict = dict(
            hgt_0 =   60., hgt_1 =  560., hgt_step=1.0,
            lat_0 =   32., lat_1 =   80., lat_step=0.50,
            lon_0 = -100., lon_1 =  -40., lon_step=0.50)

#    iono_dict = dict(
#            hgt_0 =   60., hgt_1 =  560., hgt_step=1.0,
#            lat_0 =   32., lat_1 =   80., lat_step=0.10,
#            lon_0 = -100., lon_1 =  -40., lon_step=0.10)

    # Lat/lon ctr from MUSIC calculation.
    dct = dict(
        radar                                       = 'fhe',
        lat_ctr                                     = 42.2,
        lon_ctr                                     = -95.1,
        ngates                                      = 100,
        freq                                        = 10.5,
        sTime                                       = datetime.datetime(2012,12,21,17),
        eTime                                       = datetime.datetime(2012,12,21,21),
        delta_t                                     = datetime.timedelta(minutes=2),
        iono_dict                                   = iono_dict,
        wave_list                                   = None,
        clear_output_dirs                           = False,
        run_get_ionosphere                          = True,
        load_background_ionosphere_from_pickle      = True,
##        background_ionosphere_source_for_symlink    = 'mstid_data/raytrace/fhe_20121221_paper2_14500kHz_lambda451km',
        background_ionosphere_source_for_symlink    = 'output/raytrace/fhe_20121221_background',
        load_perturbed_ionosphere_from_pickle       = False,
        save_background_ionosphere_to_pickle        = False,
        save_perturbed_ionosphere_to_pickle         = False,
        run_raytrace                                = True,
        subproc                                     = True,  # Use subprocess; this is needed to run (most) multiproc.
        multiproc                                   = True,
        nprocs                                      = 4,
        base_dir                                    = 'output/raytrace'
    )

    run_configs = []
#    run_configs.append({'freq':10.5})
    run_configs.append({'freq':14.5})

    azms_0          = np.linspace(0,360,8,endpoint=False)
    wave_defaults   = {'modulation':0.30,'rng':8000}

#    waves       = []
#    waves.append({'lambda_h':451,'azm':159.+180.,'T_minutes':40})
#    waves.append({'lambda_h':198,'azm':152.+180.,'T_minutes':40})

    beam_list_fnames    = []
    beam_list_out = 'beam_list.txt'
    if os.path.exists(beam_list_out):
        os.remove(beam_list_out)

    for run_config in run_configs:
        run_dict = dct.copy()
        run_dict.update(run_config)

        radar   = run_dict.get('radar')
        sTime   = run_dict.get('sTime')
        eTime   = run_dict.get('eTime')
        lat_ctr = run_dict.get('lat_ctr')
        lon_ctr = run_dict.get('lon_ctr')

#        for wave in waves:
        for azm_0 in azms_0:
            wave    = {'lambda_h':200,'T_minutes':40}

            this_wave   = wave_defaults.copy()
            this_wave.update(wave)

            rng         = this_wave.get('rng')
            lambda_h    = this_wave.get('lambda_h')

            wave_list = []
            src_lat, src_lon,azm    = get_source_loc(rng,azm_0,radar,sTime,lat_ctr=lat_ctr,lon_ctr=lon_ctr,azm_plus_boresite=False)
            this_wave['src_lat']    = src_lat
            this_wave['src_lon']    = src_lon
            this_wave['azm']        = azm

            wave_list.append(this_wave)

            run_dict['wave_list']   = wave_list

            txt =[]
#            txt.append('paper2')
            txt.append('azm{:03d}deg'.format(int(azm)))
            txt.append('{:.0f}kHz'.format(run_dict['freq']*1000.))
            txt.append('lambda{:.0f}km'.format(lambda_h))
            run_dict['event_note'] = '_'.join(txt)

            info_dict   = run_mstid_rt(**run_dict)
            

print "I'm done!!! :-)"
