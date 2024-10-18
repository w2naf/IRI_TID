#!/usr/bin/env python
# This is a derivative of music_and_classify.py.  No classification is done 
# here, as that requires a view of an entire MSTID season.
#
# This script will run MUSIC processing on raytrace-simulated data.  To
# generate said raytraced-simulated data, please see run_raytrace.py.

import sys
import os
import datetime
import subprocess

import numpy as np

import matplotlib
matplotlib.use('Agg')

import multiprocessing

import mstid
from mstid import run_helper

def eventize(event,base_dir='mstid_data/music_no_mongo'):
    events  = []
    defaults = {}
    defaults['db_name'] = None

    run_level_dicts = []

    tmp = {}
    tmp['data_path']            = os.path.join(base_dir,'fft')
    tmp['hanning_window_space'] = False
    tmp['process_level']        = 'fft'
    run_level_dicts.append(tmp)

#    tmp = {}
#    tmp['data_path']            = os.path.join(base_dir,'music')
#    tmp['hanning_window_space'] = True
#    tmp['process_level']        = 'music'
##    tmp['filter_cutoff_high']   = 0.00075
##    tmp['gate_limits']          = (21,35)
#    tmp['auto_range_on']        = True
#    run_level_dicts.append(tmp)

    for run_level_dict in run_level_dicts:
        event_update    = event.copy()
        event_update.update(defaults)
        event_update.update(run_level_dict)
        events.append(event_update)

    return events


# User-Defined Run Parameters Go Here. #########################################
nprocs              = 4
multiproc           = True

#******************************************************************************#
# No User Input Below This Line ***********************************************#
#******************************************************************************#
events_0 = []
#events_0.append({'radar':'fhe','sTime':datetime.datetime(2012,12,21,18),'eTime':datetime.datetime(2012,12,21,20)})

srcPaths    = []
srcPaths.append('http://sd-work1.ece.vt.edu/data/mstid/statistics/music_scripts/mstid_data/raytrace/fhe_20121221_paper2_14500kHz_lambda451km/000_general_output/fhe_20121221_paper2_14500kHz_lambda451km_beam_list.p')

for srcPath in srcPaths:
    srcPath = srcPath.split('http://sd-work1.ece.vt.edu')[-1]
    events_0.append({'radar':'fhe','sTime':datetime.datetime(2012,12,21,14),'eTime':datetime.datetime(2012,12,21,22),'srcPath':srcPath})
#    events_0.append({'radar':'fhe','sTime':datetime.datetime(2012,12,21,18),'eTime':datetime.datetime(2012,12,21,20),'srcPath':srcPath})

events = []
for event in events_0:
    event_code  = (event['srcPath'].split('/')[-1]).split('_beam_list.p')[0]
    base_dir    = os.path.join('mstid_data','music_no_mongo',event_code)
    events += eventize(event,base_dir=base_dir)

# Prepare initial_param.json files #############################################
init_files  = []
for inx,event in enumerate(events):
    prefix  = '{:03d}_'.format(inx)
    init_files.append(mstid.generate_initial_param_file(event,prefix=prefix))

if multiproc:
    if len(init_files) > 0:
        pool = multiprocessing.Pool(nprocs)
        pool.map(run_helper.run_init_file,init_files)
        pool.close()
        pool.join()
else:
    for init_file in init_files:
        cmd = ['./run_single_event.py',init_file]
        print ' '.join(cmd)
        run_helper.run_music_init_param_file(init_file)

print "I'm done!"
