#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import numpy as np
import os
#import contextlib
from astropy import units as u
import pandas as pd


def save_psd(root_dir, df_psd, df_lstar, targets, times, flag_model, time_avg, type):
    target_Ks, target_mus = targets
    start_time, end_time = times
    stime = pd.to_datetime(start_time)

    if type == 'day':
        save_dir = os.path.join(root_dir, 'Days', f'{stime.year}', f'{stime.month}', f'{stime.day}', '')
    elif type == 'storm':
            save_dir = os.path.join(root_dir, 'Storms', start_time.split()[0].replace('-', '_'), '')
    else:
        print('type not valid')
        return

    if not os.path.exists(os.path.dirname(save_dir)):
        try:
            os.makedirs(os.path.dirname(save_dir))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise

    time_array = df_psd.index
    for K in target_Ks:
        for mu in target_mus:
            psd = df_psd[str(K)][str(mu)].tolist()
            psd = u.Quantity(psd)
            lstar = df_lstar[str(K)][str(mu)].tolist()
            save_name1 = f'{flag_model}_{time_avg}_'
            save_name2 = f'K_{str(K).replace('(1/2)', '')}_mu_{str(mu).replace(' / ', '_')}.txt'.replace(' ', '_')
            save_name = save_name1 + save_name2
            save = os.path.join(save_dir, save_name)
            str_psd = 'psd'
            str_lstar = 'lstar'

            df_save = pd.DataFrame(columns = ['time', str_psd, str_lstar])
            df_save['time'] = time_array
            df_save[str_psd] = psd.value
            df_save[str_lstar] = lstar

            header_lines = [f'# Start time: {start_time}', f'# End time: {end_time}',
                            f'# K = {str(K)} ', f'# mu = {str(mu)}',
                            f'# Time avg: {time_avg}', f'# Model: {flag_model}',
                            f'# PSD units: {psd.unit}', '']
            with open(save, "w") as f:
                for line in header_lines:
                    f.write(line + "\n")
                df_save.to_csv(f, index=False)
    return df_save


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
