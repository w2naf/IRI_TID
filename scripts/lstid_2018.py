#!/usr/bin/env python
import datetime 
import os

import matplotlib
matplotlib.use('Agg')

# The GeographicLib is likely more accurate than geooack.
# Geographiclib - https://geographiclib.sourceforge.io/Python/2.0/
# mamba install conda-forge::geographiclib
from geographiclib.geodesic import Geodesic
geod = Geodesic.WGS84

import pydarn
import ionolib

output_dir  = 'output'
ionolib.gen_lib.prep_dirs({0:output_dir},clear_output_dirs=True,php=False)

profile_dir= 'output/profiles'
ionolib.gen_lib.prep_dirs({0:profile_dir},php=False)

map_dir= 'output/maps'
ionolib.gen_lib.prep_dirs({0:map_dir},php=False)

kw_args             = {}
#kw_args['engine']   = 'PyIRI'  # Victoria Forsythe's PyIRI (https://github.com/victoriyaforsythe/PyIRI)
kw_args['engine']   = 'iri2016' # Michael Hirsch's IRI2016 Python Wrapper (https://github.com/space-physics/iri2016)
kw_args['sDate']    = datetime.datetime(2018,12,15,20)
kw_args['eDate']    = datetime.datetime(2018,12,15,20)
kw_args['hgt_0']    =    0.0
kw_args['hgt_1']    =  600.0
kw_args['hgt_step'] =    3.0
kw_args['lat_0']    =   30.0
kw_args['lat_1']    =   60.0
kw_args['lon_0']    = -110.0
kw_args['lon_1']    =  -60.0
kw_args['lat_step'] =    0.10
kw_args['lon_step'] =    0.10

#kw_args['lat_0']    =   -90.
#kw_args['lat_1']    =    90.
#kw_args['lon_0']    = -180.0
#kw_args['lon_1']    =  180.0

print('Generating 3d Ionosphere...')
iono = ionolib.iono_grid.iono_3d(**kw_args)

print('Adding in TID...')
wave_list = []
#wave_list.append(dict(src_lat=40.679917,src_lon=-105.040944,amplitude=0.50,lambda_h=250,T_minutes=15))
#wave_list.append(dict(src_lat=70.,src_lon= -70.,amplitude=0.50,lambda_h=300,T_minutes=15,advance_minutes=5))
wave_list.append(dict(src_lat=60.,src_lon= 112.,amplitude=0.50,lambda_h=1000,T_minutes=120,advance_minutes=5))
iono.generate_wave(wave_list)

print('Generating ionospheric profile along chosen path...')
#radar = 'fhe'
#hdw_data = pydarn.read_hdw_file(radar,kw_args['sDate'])
#tx_lat   = hdw_data.geographic.lat
#tx_lon   = hdw_data.geographic.lon
#boresite = hdw_data.boresight.physical

radar = 'TX'
tx_lat   =  30.
tx_lon   = -85.
boresite = 0.

rx_dct   = geod.Direct(tx_lat, tx_lon, boresite, 3000e3)
rx_lat   = rx_dct['lat2']
rx_lon   = rx_dct['lon2']

prof_dct            = {}
prof_dct['tx_call'] = radar.upper()
prof_dct['tx_lat']  = tx_lat
prof_dct['tx_lon']  = tx_lon
prof_dct['rx_call'] = ''
prof_dct['rx_lat']  = rx_lat
prof_dct['rx_lon']  = rx_lon
prof_dct['range_step']  = 10.
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
import ipdb; ipdb.set_trace()
