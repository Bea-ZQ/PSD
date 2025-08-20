#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

'''
###############################################################################
### En este script están todas las funciones que me permiten crear los inputs
### para los modelos de campo magnético
################################################################################
'''

def get_omni_dates(sdate, edate):
    flag = 7 > int(sdate[5:7])
    if flag:
        sdate_omni_1h = sdate[:5] + '01-01'
    else:
        sdate_omni_1h = sdate[:5] + '07-01'
    edate_omni = pd.to_datetime(edate) + pd.Timedelta(days=1)
    sdate_omni_min = sdate[0:7] + '-01'

    print('start_date 1h: ', sdate_omni_1h)
    print('start_date min: ', sdate_omni_min)
    print('end_date: ', edate_omni)

    return sdate_omni_1h, sdate_omni_min, edate_omni


def model_variables_ect(mf_flag):
    var_ect = ['Epoch', 'Position']
    rename_ect = {'Epoch':'epoch', 'Position':'x'}
    return var_ect, rename_ect


def model_variables_omni(mf_flag, resolution):

    if mf_flag == 'T89':
        print('T89')
        var_omni_hour = ['Epoch', 'KP']
        var_omni_min = []
        rename_omni = {'Epoch':'epoch', 'KP': 'Kp'}

    elif mf_flag == 'T96':
        print('T96')
        if resolution != '1h':
            var_omni_min = ['Epoch', 'Pressure', 'BY_GSM', 'BZ_GSM']
            var_omni_hour = ['Epoch','DST']
        else:
            var_omni_hour = ['Epoch', 'Pressure', 'BY_GSM', 'BZ_GSM', 'DST']
            var_omni_min = []

        rename_omni = {'Pressure': 'Pdyn', 'BY_GSM': 'ByIMF', 'BZ_GSM': 'BzIMF',
                      'Epoch': 'epoch', 'DST': 'Dst'}

    return var_omni_hour, var_omni_min, rename_omni


def sync_omni_data(data_1h, data_min, resolution):
    omni_1h, omni_meta_1h = data_1h
    omni_min, omni_meta_min = data_min

    if resolution == '1h':
        return omni_1h, omni_meta_1h
    else:
        data_to_add = omni_1h.set_index('epoch')
        next_time = data_to_add.index[-1] + pd.Timedelta(hours=1)
        data_to_add.loc[next_time] = data_to_add.iloc[-1]
        resampled_data = data_to_add.resample(resolution).ffill()
        resampled_data.drop(index = next_time).reset_index()

        synced_omni = pd.merge(omni_min, resampled_data, on='epoch', how='left')

        synced_omni_meta = omni_meta_min | omni_meta_1h
        return synced_omni, synced_omni_meta


def time_average_ect_data(ect_info, fedu, sdate, edate, frequency):

    fedu_avg = []
    common_time = pd.date_range(sdate, pd.to_datetime(edate) + pd.Timedelta(days=1), freq = frequency)
#    print(common_time)
    bins_labels = range(len(common_time))

    ect_info['group'] = pd.cut(ect_info['epoch'], bins=common_time, labels = bins_labels[:-1])
    ect_info_avg = ect_info.groupby('group', observed='False').mean()
#    print(ect_info_avg)

    for i in range(len(ect_info_avg)):
#        print(i)
        index_bool_i = ect_info['group']==i
    #    print(index_bool_i)
        fedus_i = fedu[index_bool_i]
#        print(np.shape(fedus_i))
#        fedu_i_avg = np.nanmean(fedus_i, axis=0)
        fedu_i_avg = fedus_i.mean(axis=0)
#        print(np.shape(fedu_i_avg))
        fedu_avg.append(fedu_i_avg)
#        print('-----------')
#    print(len(fedu_avg))
#    print(type(fedu_avg))
    fedu_avg = np.array(fedu_avg)
#    print(type(fedu_avg))
#    print(np.shape(fedu_avg))
    return ect_info_avg, fedu_avg


def sync_ect_data(rept_info, mageis_info):
    combined_ect_info = pd.concat([rept_info, mageis_info], axis='index')

    averaged_ect_info = combined_ect_info.groupby('group', observed = True).mean()

    return averaged_ect_info


def filter_dates(df, sdate, edate, delta=0):
    # Filtrar el DataFrame entre las fechas especificadas

    sdate2= pd.to_datetime(sdate)
    if delta ==0:
        edate2 = pd.to_datetime(edate)
    else:
        edate2 = pd.to_datetime(edate) + delta

    mask = (df['epoch'] >= sdate2) & (df['epoch'] <= edate2)
    df_filtered = df[mask]
    #    df_filtered = df_filtered.reset_index(drop=True)
    return df_filtered


def get_GEO_coordinates(df_posit, Re):
    dates = df_posit['epoch']
    x_geo = df_posit['x1']/Re
    y_geo = df_posit['x2']/Re
    z_geo = df_posit['x3']/Re

    return dates, np.array(x_geo), np.array(y_geo), np.array(z_geo)


def dicts_x_input(dates, xs_geo, ys_geo, zs_geo):

    list_x_inputs = []
    N = len(dates)
    for i in range(N):
#        print(i)
        x_input = {'x1':np.float64(xs_geo[i]), 'x2':np.float64(ys_geo[i]),
                   'x3':np.float64(zs_geo[i]), 'dateTime':dates.iloc[i]}
        list_x_inputs.append(x_input)
    return list_x_inputs, N


def dicts_magnetic_input(dates_ect, df_omni):
    #.loc es el índice
    #.iloc es la posición
    keys_to_iterate = [col for col in df_omni.columns if col != "epoch"]
#    print(keys_to_iterate)
    list_mag_inputs = []
    start_posit_ect = 0
    for posit_omni, timestamp_omni in enumerate(df_omni['epoch'].iloc[1:]):
        end_posit_ect = dates_ect.searchsorted(timestamp_omni, side='left')
        len_i = end_posit_ect - start_posit_ect

        dict_i = {}
        for key in keys_to_iterate:
#            print(key)
            value_i = df_omni.iloc[posit_omni][key]
            dict_i[key] = float(value_i)

        mag_inputs_i = [dict_i]*len_i
        start_posit_ect = end_posit_ect
        list_mag_inputs += mag_inputs_i

    return list_mag_inputs


def create_inputs_MF(ect_info, omni_info, Re, sdate, edate):

    # datos ect
    dates, xs_geo, ys_geo, zs_geo = get_GEO_coordinates(ect_info, Re)
    list_x_inputs, N = dicts_x_input(dates, xs_geo, ys_geo, zs_geo)

    # datos omni
    delta = pd.Timedelta(days=1)
    omni_filt = filter_dates(omni_info, sdate, edate, delta)
    print(omni_filt.iloc[-1])
    list_mag_inputs = dicts_magnetic_input(dates, omni_filt)

    return (list_x_inputs, list_mag_inputs, N)


########################### Funciones de checkeo ###############################

def check_x_inputs(list_x_inputs, ect_info, Re):
    x1_flag = True
    x2_flag = True
    x3_flag = True
    for i in range(len(list_x_inputs)):
        flag1 = list_x_inputs[i]['x1'] == ect_info.iloc[i]['x1']/Re
        flag2 = list_x_inputs[i]['x2']== ect_info.iloc[i]['x2']/Re
        flag3 = list_x_inputs[i]['x3'] == ect_info.iloc[i]['x3']/Re
        x1_flag = x1_flag&flag1
        x2_flag = x2_flag&flag2
        x3_flag = x3_flag&flag3
    print('Checking x_inputs:')
    print('1', x1_flag)
    print('2', x2_flag)
    print('3', x3_flag)


def check_mag_inputs(list_mag_inputs, omni_info, ect_info, res):
    print('Checking mag_inputs:')
    keys_to_iterate = [col for col in omni_info.columns if col != "epoch"]
    for key in keys_to_iterate:
        flag_key = True
        flag_date =True
        for i in range(len(list_mag_inputs)):
            date1 = ect_info.iloc[i]['epoch']
        #    print('date1', date1)
            value1=list_mag_inputs[i][key]
            end_pos = omni_info['epoch'].searchsorted(ect_info.iloc[i]['epoch'], side='right')
            value2 = omni_info.iloc[end_pos-1][key]
            date2 = omni_info.iloc[end_pos-1]['epoch']
        #    print('date2', date2)
        #    print(value1)
        #    print(value2)
            if not (np.isnan(value1) &np.isnan(value2)):
                flag1 = value1==value2
            flag2= date1-date2< pd.Timedelta(res)
            flag_key = flag_key&flag1
            flag_date = flag_date&flag2

        print(key, flag_key)
        print('date', flag_date)
    print('----')
