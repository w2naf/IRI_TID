#!/usr/bin/env python
import os
import glob
import pandas as pd

engine      = 'iri2016'
input_dir   = os.path.join('output',engine,'spots')

spot_csvs   = glob.glob(os.path.join(input_dir,'*.spot.csv'))
spot_csvs.sort()

df_lst  = []
for spot_csv in spot_csvs:
    try:
        tmp_df = pd.read_csv(spot_csv,parse_dates=[0])
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
df  = df.sort_values('date')
import ipdb; ipdb.set_trace()
