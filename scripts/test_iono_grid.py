#!/usr/bin/env python
import datetime 

import matplotlib
matplotlib.use('Agg')

import ionolib

output_dir  = 'output'
ionolib.gen_lib.prep_dirs({0:output_dir},clear_output_dirs=True,php=False)

profile_dir= 'output/profiles'
ionolib.gen_lib.prep_dirs({0:profile_dir},php=True)

map_dir= 'output/maps'
ionolib.gen_lib.prep_dirs({0:map_dir},php=True)

kw_args             = {}
kw_args['sDate']    = datetime.datetime(2016,6,21)
kw_args['eDate']    = datetime.datetime(2016,6,21)
kw_args['hgt_0']    =   50.0
kw_args['hgt_1']    =  451.0
kw_args['hgt_step'] =    5.0
kw_args['lat_0']    =   30.0
kw_args['lat_1']    =   55.0
kw_args['lon_0']    = -110.0
kw_args['lon_1']    =  -70.0
kw_args['lat_step'] =    1.00
kw_args['lon_step'] =    1.00

iono = ionolib.iono_grid.iono_3d(**kw_args)
wave_list = []
wave_list.append(dict(src_lat=40.679917,src_lon=-105.040944,amplitude=0.50,lambda_h=250,T_minutes=15))
#wave_list.append(dict(src_lat=60.,src_lon= -80.,amplitude=0.50,lambda_h=250,T_minutes=15))
iono.generate_wave(wave_list)

prof_dct            = {}
prof_dct['tx_call'] = 'WWV'
prof_dct['tx_lat']  =   40.679917
prof_dct['tx_lon']  = -105.040944
prof_dct['rx_call'] = 'W2NAF'
prof_dct['rx_lat']  =   41.335116
prof_dct['rx_lon']  =  -75.600692
iono.generate_tx_rx_profile(**prof_dct)
iono.plot_profiles(output_dir=profile_dir)

#xlim    = (-180,180)
#ylim    = (-90,90)

xlim    = (-130,-56)
ylim    = (20,55)
iono.plot_maps(output_dir=map_dir,xlim=xlim,ylim=ylim)

import ipdb; ipdb.set_trace()
