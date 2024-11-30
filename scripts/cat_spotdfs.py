#!/usr/bin/env python
import os
import glob
import datetime
import tqdm

import pandas as pd
import numpy as np

engine      = 'iri2016'
base_dir    = os.path.join('output',engine)
input_dir   = os.path.join(base_dir,'spots')

spot_csvs   = glob.glob(os.path.join(input_dir,'*.spot.csv'))
spot_csvs.sort()

df_lst  = []
for spot_csv in spot_csvs:
    try:
#        tmp_df = pd.read_csv(spot_csv,parse_dates=[0])
        tmp_df = pd.read_csv(spot_csv)
        df_lst.append(tmp_df)
        print(f'LOADED: {spot_csv}')
    except:
        print(f'SKIPPING: {spot_csv}')
        pass

n_loaded    = len(df_lst)
n_csvs      = len(spot_csvs)
pct_loaded  = n_loaded/n_csvs*100
print()
print(f'Loaded {n_loaded}/{n_csvs} ({pct_loaded:.0f}%) Spot CSVs')
df  = pd.concat(df_lst,ignore_index=True)

# For some inexplicable reason, using pd.read_csv(spot_csv,parse_dates[0]) does not work reliably.
# However, looping through and doing each pd.Timestamp conversion separately as below does.
for rinx, row in tqdm.tqdm(df.iterrows(),total=len(df),desc='Parsing Dates',dynamic_ncols=True):
    ts = pd.Timestamp(row['date'])
    df.loc[rinx,'date'] = ts

df  = df.sort_values('date')

sDate   = df['date'].min()
eDate   = df['date'].max()
sDate_str   = sDate.strftime('%Y%m%d.%H%M')
eDate_str   = eDate.strftime('%Y%m%d.%H%M')
out_fname   = f'{sDate_str}-{eDate_str}_{engine}.spot.csv'
out_fpath   = os.path.join(base_dir,out_fname)
df.to_csv(out_fpath,index=False)

import ipdb; ipdb.set_trace()
