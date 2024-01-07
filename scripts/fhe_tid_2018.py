#!/usr/bin/env python
import datetime 

import matplotlib
matplotlib.use('Agg')

import pydarn

import ionolib

output_dir  = 'output'
ionolib.gen_lib.prep_dirs({0:output_dir},clear_output_dirs=True,php=False)

profile_dir= 'output/profiles'
ionolib.gen_lib.prep_dirs({0:profile_dir},php=True)

map_dir= 'output/maps'
ionolib.gen_lib.prep_dirs({0:map_dir},php=True)

kw_args             = {}
kw_args['sDate']    = datetime.datetime(2018,12,10,18,30)
kw_args['eDate']    = datetime.datetime(2018,12,10,18,30)
kw_args['hgt_0']    =   50.0
kw_args['hgt_1']    =  451.0
kw_args['hgt_step'] =    5.0
kw_args['lat_0']    =   30.0
kw_args['lat_1']    =   55.0
kw_args['lon_0']    = -110.0
kw_args['lon_1']    =  -70.0
kw_args['lat_step'] =    1.00
kw_args['lon_step'] =    1.00

#kw_args['lat_0']    =   -90.
#kw_args['lat_1']    =    90.
#kw_args['lon_0']    = -180.0
#kw_args['lon_1']    =  180.0

iono = ionolib.iono_grid.iono_3d(**kw_args)
wave_list = []
#wave_list.append(dict(src_lat=40.679917,src_lon=-105.040944,amplitude=0.50,lambda_h=250,T_minutes=15))
wave_list.append(dict(src_lat=50.,src_lon= -60.,amplitude=0.50,lambda_h=250,T_minutes=15))
#iono.generate_wave(wave_list)

radar = 'fhe'
hdw_data = pydarn.read_hdw_file(radar,kw_args['sDate'])
tx_lat   = hdw_data.geographic.lat
tx_lon   = hdw_data.geographic.lon
boresite = hdw_data.boresight.physical


rx_dct   = ionolib.geopack.calcDistPnt(tx_lat,tx_lon,0.,az=boresite,el=0,dist=2000.)
rx_lat   = rx_dct['distLat']
rx_lon   = rx_dct['distLon']

prof_dct            = {}
prof_dct['tx_call'] = radar.upper()
prof_dct['tx_lat']  = tx_lat
prof_dct['tx_lon']  = tx_lon
prof_dct['rx_call'] = ''
prof_dct['rx_lat']  = rx_lat
prof_dct['rx_lon']  = rx_lon
iono.generate_tx_rx_profile(**prof_dct)
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

iono.plot_maps(output_dir=map_dir,xlim=xlim,ylim=ylim)
#iono.plot_maps_ortho(output_dir=map_dir,xlim=xlim,ylim=ylim)

import ipdb; ipdb.set_trace()
