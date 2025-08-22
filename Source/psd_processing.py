#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import numpy as np
import os
#import contextlib
from astropy import units as u
import pandas as pd



def read_psd(root_dir, K, mu, start_time, flag_model, time_avg, type):
    stime = pd.to_datetime(start_time)
    if type == 'day':
        data_dir = os.path.join(root_dir, 'Days', f'{stime.year}', f'{stime.month}', f'{stime.day}', '')
    elif type == 'storm':
        #ESTO NSE PUEDE HACER COMO EL CASO DE ARRIBA
        data_dir = os.path.join(root_dir, 'Storms', f'{stime.year}_{stime.month}_{stime.day}', '')
    else:
        print('type not valid')
        return

    name_1 = f'{flag_model}_{time_avg}_'
    name_2 = f'K_{str(K).replace('(1/2)', '')}_mu_{str(mu).replace(' / ', '_')}.txt'.replace(' ', '_')
    name = name_1 + name_2
    file = os.path.join(data_dir, name)

    with open(file, 'r') as f:
        for line in f:
            if 'PSD units' in line:
                unit_psd = line.split(':')[1].strip()
                break

    df = pd.read_csv(file, header = 7)
    df['time'] = pd.to_datetime(df['time'])
    return df, unit_psd


def clean_psd(df):
    df_clean = df.dropna(axis = 0, how='any')
    df_clean = df_clean[df_clean['psd'] >= 10e-13]
    df_clean = df_clean[df_clean['psd'] <= 10]
    return df_clean
