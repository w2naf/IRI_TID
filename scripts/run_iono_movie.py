#!/usr/bin/env python
import os
import datetime
import multiprocessing
import subprocess

def gen_cmd(rd):
    """
    Generates commands in the form of:
    cmd = './iono_frame.py --time="2018-12-15 12:00" --engine=PyIRI'
    """

    cmd = ['./iono_frame.py']
    for opt, val in rd.items():
        if opt in ['time','time_0']:
            val = f'"{val}"'

        txt = f'--{opt}={val}'
        cmd.append(txt)
    cmd = ' '.join(cmd)
    return cmd

def run_cmd(cmd):
    return subprocess.call(cmd, shell=True)

multiproc   = True
engine      = 'PyIRI'
sTime       = datetime.datetime(2018,12,15,12)
eTime       = sTime + datetime.timedelta(hours=12)
dt          = datetime.timedelta(minutes=10)

times   = [sTime]
while times[-1] < eTime:
    times.append(times[-1] + dt)

run_dicts = []
for time in times:
    rd = {}
    rd['time']      = time
    rd['time_0']    = times[0]
    rd['engine']    = engine
    run_dicts.append(rd)

cmds = [gen_cmd(rd) for rd in run_dicts]

if multiproc:
    count   = multiprocessing.cpu_count()
    pool    = multiprocessing.Pool(processes=count)
    pool.map(run_cmd,cmds)
else:
    for cmd in cmds:
        print(cmd)
        run_cmd(cmd)

