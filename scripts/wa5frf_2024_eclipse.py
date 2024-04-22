#!/usr/bin/env python
import datetime 
import os
import pandas as pd

import matplotlib
matplotlib.use('Agg')

# The GeographicLib is likely more accurate than geooack.
# Geographiclib - https://geographiclib.sourceforge.io/Python/2.0/
# mamba install conda-forge::geographiclib
from geographiclib.geodesic import Geodesic
geod = Geodesic.WGS84

import ionolib

import plotLAP

output_dir  = 'eclipse2024'
ionolib.gen_lib.prep_dirs({0:output_dir},clear_output_dirs=True,php=False)

profile_dir = os.path.join(output_dir,'profiles')
ionolib.gen_lib.prep_dirs({0:profile_dir},php=False)

map_dir     = os.path.join(output_dir,'maps')
ionolib.gen_lib.prep_dirs({0:map_dir},php=False)

raytrace_dir     = os.path.join(output_dir,'raytraces')
ionolib.gen_lib.prep_dirs({0:raytrace_dir},php=False)

kw_args             = {}
#kw_args['engine']   = 'PyIRI'  # Victoria Forsythe's PyIRI (https://github.com/victoriyaforsythe/PyIRI)
kw_args['engine']   = 'iri2016' # Michael Hirsch's IRI2016 Python Wrapper (https://github.com/space-physics/iri2016)
kw_args['sDate']    = datetime.datetime(2020,4,8,19,30) # iri2016 does not seem to give results after 2020.
kw_args['eDate']    = datetime.datetime(2020,4,8,19,30)
kw_args['hgt_0']    =    0.0
kw_args['hgt_1']    =  600.0
kw_args['hgt_step'] =    3.0
kw_args['lat_0']    =   25.0
kw_args['lat_1']    =   50.0
kw_args['lon_0']    = -120.0
kw_args['lon_1']    =  -70.0
kw_args['lat_step'] =    0.50
kw_args['lon_step'] =    0.50

#kw_args['lat_0']    =   -90.
#kw_args['lat_1']    =    90.
#kw_args['lon_0']    = -180.0
#kw_args['lon_1']    =  180.0

print('Generating 3d Ionosphere...')
iono = ionolib.iono_grid.iono_3d(**kw_args)

#print('Adding in TID...')
#wave_list = []
##wave_list.append(dict(src_lat=40.679917,src_lon=-105.040944,amplitude=0.50,lambda_h=250,T_minutes=15))
#wave_list.append(dict(src_lat=70.,src_lon= -70.,amplitude=0.50,lambda_h=300,T_minutes=15,advance_minutes=5))
#iono.generate_wave(wave_list)

print('Generating ionospheric profile along chosen path...')

prof_dct            = {}
prof_dct['tx_call'] = 'WA5FRF'
prof_dct['tx_lat']  =  29.574936802385718
prof_dct['tx_lon']  = -98.8871154108325
prof_dct['rx_call'] = 'W3USR'
prof_dct['rx_lat']  =  41.40528906979763
prof_dct['rx_lon']  = -75.65787801145272 
prof_dct['range_step']  = 10.
prof_dct['max_range']   = 3000. + prof_dct['range_step']
iono.generate_tx_rx_profile(**prof_dct)

print('Saving ionospheric profile to netcdf in {!s}'.format(profile_dir))
iono.profiles_to_netcdf(output_dir=profile_dir)

print('Plotting ionospheric profile to PNGs in {!s}'.format(profile_dir))
iono.plot_profiles(output_dir=profile_dir)

## World
#xlim    = (-180,180)
#ylim    = (-90,90)

# CONUS
xlim    = (-130,-56)
ylim    = (20,55)

## CONUS + Canada
#xlim    = (-130,-56)
#ylim    = (20,80)

print('Saving 3D ionospheric dataset to netcdf in {!s}'.format(map_dir))
iono.iri_dataset.to_netcdf(os.path.join(map_dir,iono.fname))

print('Plotting ionospheric map to PNGs in {!s}'.format(map_dir))
iono.plot_maps(output_dir=map_dir,xlim=xlim,ylim=ylim)
#iono.plot_maps_ortho(output_dir=map_dir,xlim=xlim,ylim=ylim)

# Raytrace profiles
freq            = 14.300
plot_end_range  = None
for profile_key,profile in iono.profiles.items():
    UT          = pd.to_datetime(profile['date'].values[0])
    date_str    = UT.strftime('%Y%m%d_%H%MUT')
    png_fname   = '{!s}_{!s}_{:0.0f}kHz_raytrace.png'.format(date_str,profile.attrs['fname_base'],freq*1000)
    png_fpath   = os.path.join(raytrace_dir,png_fname)
    plotLAP.main(profile,freq=freq,png_fpath=png_fpath,plot_end_range=plot_end_range)
    import ipdb; ipdb.set_trace()
