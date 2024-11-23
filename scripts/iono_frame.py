#!/usr/bin/env python
# This script is used to generate the LSTID-perturbed IRI profile used in the 
# Figure 2d raytrace for Frissell et al. (2024) GRL, submitted.
# Nathaniel A. Frissell - 6 September 2024

usage = """
    Usage:

    iono_frame.py --time=<YYYY-MM-DDTHH:mm>     \
                  --time_0=<YYYY-MM-DDTHH:mm>   \
                  --engine=<PyIRI or iri2016>

    where:
    
    --time=<YYYY-MM-DDTHH:mm>   - Datetime for ionosphere run

    --time_0=<YYYY-MM-DDTHH:mm> - Zero reference for phase advance of added TIDs

    --engine=<PyIRI or iri2016> - IRI engine. Options are:
        'PyIRI':    Victoria Forsythe's PyIRI (https://github.com/victoriyaforsythe/PyIRI)
        'iri2016':  Michael Hirsch's IRI2016 Python Wrapper (https://github.com/space-physics/iri2016)
    """

import sys
import os
import getopt
import datetime 

import matplotlib
matplotlib.use('Agg')
import pandas as pd

# The GeographicLib is likely more accurate than geooack.
# Geographiclib - https://geographiclib.sourceforge.io/Python/2.0/
# mamba install conda-forge::geographiclib
from geographiclib.geodesic import Geodesic
geod = Geodesic.WGS84

import ionolib

# Parse Command Line
arglist = ''
longarglist = ['time=',
               'time_0=',
               'engine=']

optlist, args = getopt.getopt(sys.argv[1:], arglist, longarglist)

# Set Default Values
time    = datetime.datetime(2018,12,15,12)
time_0  = None
#engine  = 'PyIRI'  # Victoria Forsythe's PyIRI (https://github.com/victoriyaforsythe/PyIRI)
engine  = 'iri2016' # Michael Hirsch's IRI2016 Python Wrapper (https://github.com/space-physics/iri2016)

for opt in optlist:
    if opt[0] == '--time':
        time    = pd.to_datetime(opt[1])
    elif opt[0] == '--time_0':
        time_0  = pd.to_datetime(opt[1])
    elif opt[0] == '--engine':
        engine = opt[1]
    else:
        raise ValueError('Illegal option %s\n%s' % (opt[0], usage))

# verify that no regular arguments were passed in
if len(args) != 0:
    raise ValueError('This command does not accept any arguments without options - may be due to illegal spaces in the command')

if time_0 is None:
    time_0 = time

output_dir  = os.path.join('output',engine)
profile_dir = os.path.join(output_dir,'profiles')
ionolib.gen_lib.prep_dirs({0:profile_dir},php=False)

map_dir     = os.path.join(output_dir,'maps')
ionolib.gen_lib.prep_dirs({0:map_dir},php=False)

# Create a log file.
log_dir     = os.path.join(output_dir,'logs')
ionolib.gen_lib.prep_dirs({0:log_dir},php=False)
log_fname   = []
log_fname.append(time.strftime('%Y%m%dT%H%Mz'))
log_fname.append(time_0.strftime('T0-%Y%m%dT%H%Mz'))
log_fname.append(engine)
log_fname.append('iono_frame')
log_fname   = '_'.join(log_fname)+'.log'
log_fpath   = os.path.join(log_dir,log_fname)
print(f'Logging to {log_fpath}')
with open(log_fpath,'w') as log_fl:
    log_fl.write(f'Log File: {log_fname}\n')

kw_args             = {}
kw_args['engine']   = engine
kw_args['sDate']    = time
kw_args['eDate']    = time
kw_args['hgt_0']    =    0.0
kw_args['hgt_1']    =  600.0
kw_args['hgt_step'] =    3.0

kw_args['lat_0']    =   20.0
kw_args['lat_1']    =   80.0
kw_args['lon_0']    = -130.0
kw_args['lon_1']    =  -50.0
kw_args['lat_step'] =    0.50
kw_args['lon_step'] =    0.50

#kw_args['lat_0']    =   -90.
#kw_args['lat_1']    =    90.
#kw_args['lon_0']    = -180.0
#kw_args['lon_1']    =  180.0

print('Generating 3d Ionosphere...')
iono = ionolib.iono_grid.iono_3d(**kw_args)

print('Adding in TID...')
advance_minutes = (time - time_0).total_seconds()/60.
#wave_list = []
#wave_list.append(dict(src_lat=40.679917,src_lon=-105.040944,amplitude=0.50,lambda_h=250,T_minutes=15))
#wave_list.append(dict(src_lat=70.,src_lon= -70.,amplitude=0.50,lambda_h=300,T_minutes=15,advance_minutes=5))
#wave_list.append(dict(src_lat=60.,src_lon= 112.,amplitude=0.50,lambda_h=1000,T_minutes=120,advance_minutes=advance_minutes))

wave_dct    = dict(src_lat=60.,src_lon= 112.,amplitude=0.50,lambda_h=1000,T_minutes=120,advance_minutes=advance_minutes)
iono.generate_wave([wave_dct])

print('Generating ionospheric profile along chosen path...')

paths_fname = '20181215_14000-14350kHz_montePaths.csv'
df_paths    = pd.read_csv(paths_fname,comment='#')


prof_dcts   = []
for rinx, row in df_paths.iterrows():
    prof_dct                = {}
    prof_dct['tx_call']     = row['call_sign_tx']
    prof_dct['tx_lat']      = row['txlat']
    prof_dct['tx_lon']      = row['txlon']
    prof_dct['rx_call']     = row['call_sign_rx']
    prof_dct['rx_lat']      = row['rxlat']
    prof_dct['rx_lon']      = row['rxlon']
    prof_dct['range_step']  = 10.
    prof_dct['max_range']   = 3000.
    prof_dct['interp_type'] = 'nearest'

    attrs                   = {}
    attrs['engine']         = engine
    attrs['tx_rx_pthlen']   = row['pthlen']
    attrs['tx_rx_latcen']   = row['latcen']
    attrs['tx_rx_loncen']   = row['loncen']
    attrs['tx_rx_azm']      = row['azm']
    attrs['message']        = f'Wave: {wave_dct}'

    prof_dct['attrs']       = attrs
    prof_dcts.append(prof_dct)

for prof_dct in prof_dcts:
    print(f'{datetime.datetime.now()}: iono.generate_tx_rx_profile({prof_dct})')
    with open(log_fpath,'a') as log_fl:
        log_fl.write(f'{datetime.datetime.now()}: iono.generate_tx_rx_profile({prof_dct})\n')

    iono.generate_tx_rx_profile(**prof_dct)

print('Saving ionospheric profile to netcdf in {!s}'.format(profile_dir))
iono.profiles_to_netcdf(output_dir=profile_dir)

print('Plotting ionospheric profile to PNGs in {!s}'.format(profile_dir))
iono.plot_profiles(output_dir=profile_dir)

## World
#xlim    = (-180,180)
#ylim    = (-90,90)

## CONUS
#xlim    = (-130,-56)
#ylim    = (20,55)

# CONUS + Canada
xlim    = (-130,-56)
ylim    = (20,80)

print('Saving 3D ionospheric dataset to netcdf in {!s}'.format(map_dir))
iono.iri_dataset.to_netcdf(os.path.join(map_dir,iono.fname))

print('Plotting ionospheric map to PNGs in {!s}'.format(map_dir))
iono.plot_maps(output_dir=map_dir,xlim=xlim,ylim=ylim)
#iono.plot_maps_ortho(output_dir=map_dir,xlim=xlim,ylim=ylim)
