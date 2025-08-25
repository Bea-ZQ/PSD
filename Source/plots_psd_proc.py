#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib
#import numpy as np
#from astropy import units as u
#import os
#import contextlib
import matplotlib.dates as mdates
import pandas as pd

#matplotlib.rcParams['text.usetex'] = True
#plt.rcParams['font.family'] = 'serif'

'''
###############################################################################
### En este script están todas las funciones para visualizar que los cálculos
### que involucran invariantes adiabáticas sean correctos
################################################################################
'''

def errors_channel(ax, x, df_rmse, df_r2, channel):
    rmse= df_rmse[channel]
    r2 = df_r2[channel]

    ax.plot(x, rmse, '-o', markersize= 4, color = 'k', label = 'rmse')
    ax.plot(x, r2, '-*', markersize= 4, color = 'r', label = 'r2')


def psd_lstar(df, df_rmse, df_r2, K, mu, unit_psd, flag_error, list_channels=''):
    lstar = df['lstar'].tolist()
    psd = df['psd'].tolist()

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)
    ax.plot(lstar, psd, '-o', markersize= 4, color = 'b')

    if flag_error:
        for channel in list_channels:
            fig2, ax2 = plt.subplots()
            fig2.set_size_inches(8, 7)
            errors_channel(ax2, lstar, df_rmse, df_r2, channel)
            ax2.set_title('Errors\n'+channel, fontsize = 20)
            ax2.set_xlabel('$L^*$', fontsize=17)
            ax2.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)

    title = 'K = ' + str(K) + ',    $\\mu$ = ' + str(mu)
    ax.set_yscale('log')
    ax.set_ylabel(f'PSD [{unit_psd}]', fontsize=17)
    ax.set_xlabel('$L^*$', fontsize=17)
#    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
    ax.tick_params(axis='both', labelsize=16)
    ax.set_title(title, fontsize = 20)

    plt.show()


def psd_time(df, df_rmse, df_r2, K, mu, unit_psd, flag_error, list_channels=''):
    time = df['time'].tolist()
    psd = df['psd'].tolist()

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)
    ax.plot(time, psd, '-o', markersize= 4, color = 'b')

    if flag_error:
        for channel in list_channels:
            fig2, ax2 = plt.subplots()
            fig2.set_size_inches(8, 7)
            errors_channel(ax2, time, df_rmse, df_r2, channel)
            ax2.set_title('Errors\n'+channel, fontsize = 20)
            ax2.set_xlabel('Time', fontsize=17)
            ax2.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
            ax2.xaxis.set_major_locator(mdates.HourLocator(interval = 1))  # Localizador de días
            ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))  # Formato de fecha


    title = 'K = ' + str(K) + ',    $\\mu$ = ' + str(mu)
    ax.set_yscale('log')
    ax.set_ylabel(f'PSD [{unit_psd}]', fontsize=17)
    ax.set_xlabel('Time', fontsize=17)
    ax.tick_params(axis='both', labelsize=16)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval = 1))  # Localizador de días
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))  # Formato de fecha
    ax.set_title(title, fontsize = 20)
    plt.show()


def psd_lstar_list_inv(list_data, fixed_inv, colors, unit_psd):
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)

    for i, elem in enumerate(list_data):
        df, df_rmse, df_r2, inv = elem
        lstar = df['lstar'].tolist()
        psd = df['psd'].tolist()
        ax.plot(lstar, psd, '-o', markersize= 4, color = colors[i], label = inv)

    title = fixed_inv
    ax.set_yscale('log')
    ax.set_ylabel(f'PSD [{unit_psd}]', fontsize=17)
    ax.set_xlabel('$L^*$', fontsize=17)
    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
    ax.tick_params(axis='both', labelsize=16)
    ax.set_title(title, fontsize = 20)
    plt.show()


def psd_lstar_list_times(list_data, K, mu, colors, unit_psd):
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)
    for i, elem in enumerate(list_data):
        df, df_rmse, df_r2, time = elem
        lstar = df['lstar'].tolist()
        psd = df['psd'].tolist()
        ax.plot(lstar, psd, '-o', markersize= 4, color = colors[i],
                label = f'{time} - ')

    title = f'$K$ = {str(K)}, $\\mu$ = {str(mu)}'
    ax.set_yscale('log')
    ax.set_ylabel(f'PSD [{unit_psd}]', fontsize=17)
    ax.set_xlabel('$L^*$', fontsize=17)
    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
    ax.tick_params(axis='both', labelsize=16)
    ax.set_title(title, fontsize = 20)
    plt.show()


def psd_time_list_inv(list_data, fixed_inv, colors, unit_psd):
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)

    for i, elem in enumerate(list_data):
        df, df_rmse, df_r2, inv = elem
        time = df['time'].tolist()
        psd = df['psd'].tolist()
        ax.plot(time, psd, '-o', markersize= 4, color = colors[i], label = inv)

    title = fixed_inv
    ax.set_yscale('log')
    ax.set_ylabel(f'PSD [{unit_psd}]', fontsize=17)
    ax.set_xlabel('Time', fontsize=17)
    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
    ax.tick_params(axis='both', labelsize=16)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval = 1))  # Localizador de días
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))  # Formato de fecha
    ax.set_title(title, fontsize = 20)
    plt.show()
