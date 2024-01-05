#!/usr/bin/env python
import datetime 

import matplotlib
matplotlib.use('Agg')

import ionolib

output_dir  = 'output'
ionolib.gen_lib.prep_dirs({0:output_dir},clear_output_dirs=True,php=False)

sTime       = datetime.datetime(2012,12,25)
eTime       = datetime.datetime(2012,12,25,12)

kw_args          = {}
kw_args['hgt_0']     = 0.
kw_args['hgt_1']     = 601.
kw_args['hgt_step']  = 3
kw_args['lat_0']     = 20.
kw_args['lon_0']     = -130.
kw_args['lat_1']     = 55.
kw_args['lon_1']     = -56.
kw_args['lat_step']  = 1.0
kw_args['lon_step']  = 1.0

while sTime < eTime:
    iono = ionolib.iono_slice.iono_3d(sTime,**kw_args)

    tx_lat, tx_lon  = (44.838,-123.603)
    rx_lat, rx_lon  = (33.003,-79.420)
    iono.generate_tx_rx_profile(tx_lat,tx_lon,rx_lat,rx_lon)
    iono.generate_wave()
    iono.plot_profiles()

    iono.plot_maps()


    import ipdb; ipdb.set_trace()
    sTime += datetime.timedelta(minutes=60.)

import ipdb; ipdb.set_trace()
