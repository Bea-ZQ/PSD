#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import numpy as np
import os
#import contextlib
from astropy import units as u
import pandas as pd


def get_save_dir(root_dir, start_time, type, create=True):

    stime = pd.to_datetime(start_time)

    if type == 'day':
        save_dir = os.path.join(root_dir, 'Days', f'{stime.year}', f'{stime.month}', f'{stime.day}', '')
    elif type == 'storm':
        save_dir = os.path.join(root_dir, 'Storms', start_time.split()[0].replace('-', '_'), '')
    else:
        print('type not valid')
        return 0
    if create:
        if not os.path.exists(os.path.dirname(save_dir)):
            try:
                os.makedirs(os.path.dirname(save_dir))
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise
    return save_dir


def get_save_name(flag_model, time_avg, sat, K, mu):

    save_name1 = f'{flag_model}_{time_avg}_'
    save_name2 = f'K_{str(K).replace('(1/2)', '')}_mu_{str(mu).replace(' / ', '_')}'.replace(' ', '_')
    save_name = f'rbsp{sat}_'+ save_name1 + save_name2
    return save_name


def save_df(save, df_save, lines):
    with open(save, "w") as f:
        for line in lines:
            f.write(line + "\n")
        f.write("\n")
        df_save.to_csv(f, index=False)


def save_psd_data(root_dir, df_psd, df_lstar, df_rmse, df_r2, targets, times, flag_model, time_avg, sat, type):
    target_Ks, target_mus = targets
    start_time, end_time = times
    save_dir = get_save_dir(root_dir, start_time, type)

    if save_dir == 0:
        return

    time_array = df_psd.index
    df_rmse['time'] = df_rmse.index
    df_r2['time'] = df_r2.index
    for K in target_Ks:
        for mu in target_mus:
            psd = df_psd[str(K)][str(mu)].tolist()
            psd = u.Quantity(psd)
            lstar = df_lstar[str(K)][str(mu)].tolist()

            save_name = get_save_name(flag_model, time_avg, sat, K, mu)

            save = os.path.join(save_dir, save_name+'.txt')
            save_rmse = os.path.join(save_dir, save_name+ '_rmse.txt')
            save_r2 = os.path.join(save_dir, save_name + '_r2.txt')

            str_psd = 'psd'
            str_lstar = 'lstar'
            df_save = pd.DataFrame(columns = ['time', str_psd, str_lstar])
            df_save['time'] = time_array
            df_save[str_psd] = psd.value
            df_save[str_lstar] = lstar

            header_lines = [f'# Start time: {start_time}', f'# End time: {end_time}',
                            f'# K = {str(K)} ', f'# mu = {str(mu)}',
                            f'# Time avg: {time_avg}', f'# Model: {flag_model}',
                            f'# PSD units: {psd.unit}']

            save_df(save, df_save, header_lines)
            save_df(save_rmse, df_rmse, header_lines[:-1])
            save_df(save_r2, df_r2, header_lines[:-1])
    return


def read_psd(root_dir, K, mu, start_time, flag_model, time_avg, sat, type):
    stime = pd.to_datetime(start_time)
    data_dir = get_save_dir(root_dir, start_time, type, create=False)
    if data_dir == 0:
        return
    name = get_save_name(flag_model, time_avg, sat, K, mu)

    file = os.path.join(data_dir, name+'.txt')
    file_rmse = os.path.join(data_dir, name+'_rmse.txt')
    file_r2 = os.path.join(data_dir, name+'_r2.txt')

    with open(file, 'r') as f:
        for line in f:
            if 'PSD units' in line:
                unit_psd = line.split(':')[1].strip()
                break

    df = pd.read_csv(file, header = 7)
    df_rmse = pd.read_csv(file_rmse, skiprows=6)
    df_r2 = pd.read_csv(file_r2, skiprows=6)
    df['time'] = pd.to_datetime(df['time'])
    df_rmse['time'] = pd.to_datetime(df_rmse['time'])
    df_r2['time'] = pd.to_datetime(df_r2['time'])
    return df, unit_psd, df_rmse, df_r2


def clean_psd(df_psd, df_rmse, df_r2):
    df_psd_clean = df_psd.dropna(axis = 0, how='any')
    df_psd_clean = df_psd_clean[df_psd_clean['psd'] >= 10e-13]
    df_psd_clean = df_psd_clean[df_psd_clean['psd'] <= 10]

    df_rmse_clean = df_rmse.loc[df_psd_clean.index]
    df_r2_clean = df_r2.loc[df_psd_clean.index]

    return df_psd_clean, df_rmse_clean, df_r2_clean
