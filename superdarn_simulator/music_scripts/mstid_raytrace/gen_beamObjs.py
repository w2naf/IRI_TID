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

import pydarn
from models import raydarn
import iono_lib
from iono_lib import iono_3d

def df_log_power(rt_run,fov):
    """Compute the log power in an rt_run dataframe."""
    fit_df = rt_run['fit_df']

    # Find the dates and make sure they are actually datetime.datetime objects.
    dates   = fit_df['date'].unique()
    dates   = [pd.to_datetime(x).to_datetime() for x in dates]
    
    import ipdb; ipdb.set_trace()

def rt_run_to_beam_objects(rt_run):
    #Beam record FROM: 2011-01-01 01:00:24.789000
    #bmnum  = 7 
    #fPtr  = object 
    #fitex  = None 
    #fit  = object 
    #prm  = object 
    #recordDict  = object 
    #stid  = 33 
    #lmfit  = None 
    #exflg  = None 
    #iqflg  = None 
    #offset  = 5075345 
    #rawacf  = object 
    #lmflg  = None 
    #rawflg  = None 
    #fType  = fitacf 
    #time  = 2011-01-01 01:00:24.789000 
    #acflg  = None 
    #cp  = 153 
    #iqdat  = object 
    #fitacf  = None 
    #channel  = 1 

    network = pydarn.radar.network()
    for inx,row in rt_run_data['fit_df'].iterrows():
        myBeam = pydarn.sdio.beamData()

        myBeam.bmnum        = row['beam']   # Beam number
        myBeam.fPtr         = None          # File pointer (not used)
        myBeam.fitex        = None          # Fitex flag
    #    myBeam.fit          = object 
    #    myBeam.prm          = object 
        myBeam.recordDict   = None
        myBeam.stid         = network.getRadarByCode(row['radar'])
        myBeam.lmfit        = None 
        myBeam.exflg        = None 
        myBeam.iqflg        = None 
        myBeam.offset       = None
    #    myBeam.rawacf       = object 
        myBeam.lmflg        = None 
        myBeam.rawflg       = None 
        myBeam.fType        = 'fitacf'
        myBeam.time         = row['date']
        myBeam.acflg        = None 
        myBeam.cp           = -7372
    #    myBeam.iqdat        = object 
        myBeam.fitacf       = None 
        myBeam.channel      = 1 
    import ipdb; ipdb.set_trace()

if __name__ == "__main__":
    filename    = os.path.join('data','event','rt_run_data.p')
    with open(filename,'rb') as fl:
        rt_run_data = pickle.load(fl)

    df_log_power(rt_run_data)
    rt_run_to_beam_objects(rt_run_data)
