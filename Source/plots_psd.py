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

def psd_lstar(df, K, mu, unit_psd):
    lstar = df['lstar'].tolist()
    psd = df['psd'].tolist()

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)

    ax.plot(lstar, psd, '-o', markersize= 4, color = 'b')

    title = 'K = ' + str(K) + ',    $\\mu$ = ' + str(mu)
    ax.set_yscale('log')
    ax.set_ylabel(f'PSD [{unit_psd}]', fontsize=17)
    ax.set_xlabel('$L^*$', fontsize=17)
#    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
    ax.tick_params(axis='both', labelsize=16)
    ax.set_title(title, fontsize = 20)
    plt.show()


def psd_lstar_Ks(list_dfs, list_Ks, mu, colors, unit_psd):

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)
    for i, df in enumerate(list_dfs):
        lstar = df['lstar'].tolist()
        psd = df['psd'].tolist()
        ax.plot(lstar, psd, '-o', markersize= 4, color = colors[i], label = f'$K$ = {list_Ks[i]}')

    title = '$\\mu$ = ' + str(mu)
    ax.set_yscale('log')
    ax.set_ylabel(f'PSD [{unit_psd}]', fontsize=17)
    ax.set_xlabel('$L^*$', fontsize=17)
    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
    ax.tick_params(axis='both', labelsize=16)
    ax.set_title(title, fontsize = 20)
    plt.show()


def psd_lstar_mus(list_dfs, list_mus, K, colors, unit_psd):

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)
    for i, df in enumerate(list_dfs):
        lstar = df['lstar'].tolist()
        psd = df['psd'].tolist()
        ax.plot(lstar, psd, '-o', markersize= 4, color = colors[i], label = f'$\\mu$ = {list_mus[i]}')

    title = f'$K$ = {str(K)}'
    ax.set_yscale('log')
    ax.set_ylabel(f'PSD [{unit_psd}]', fontsize=17)
    ax.set_xlabel('$L^*$', fontsize=17)
    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
    ax.tick_params(axis='both', labelsize=16)
    ax.set_title(title, fontsize = 20)
    plt.show()


def psd_lstar_times(list_dfs, list_times, mu, K, colors, unit_psd):
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)
    for i, df in enumerate(list_dfs):
        lstar = df['lstar'].tolist()
        psd = df['psd'].tolist()
        ax.plot(lstar, psd, '-o', markersize= 4, color = colors[i],
                label = f'{list_times[i]} - ')

    title = f'$K$ = {str(K)}, $\\mu$ = {str(mu)}'
    ax.set_yscale('log')
    ax.set_ylabel(f'PSD [{unit_psd}]', fontsize=17)
    ax.set_xlabel('$L^*$', fontsize=17)
    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
    ax.tick_params(axis='both', labelsize=16)
    ax.set_title(title, fontsize = 20)
    plt.show()


def psd_time(df, K, mu, unit_psd):
    time = df['time'].tolist()
    psd = df['psd'].tolist()
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)

    ax.plot(time, psd, '-o', markersize= 4, color = 'b')

    title = 'K = ' + str(K) + ',    $\\mu$ = ' + str(mu)
    ax.set_yscale('log')
    ax.set_ylabel(f'PSD [{unit_psd}]', fontsize=17)
    ax.set_xlabel('Time', fontsize=17)
#    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
    ax.tick_params(axis='both', labelsize=16)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval = 1))  # Localizador de días
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))  # Formato de fecha
    ax.set_title(title, fontsize = 20)
    plt.show()


def psd_time_Ks(list_dfs, list_Ks, mu, colors, unit_psd):

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)
    for i, df in enumerate(list_dfs):
        time = df['time'].tolist()
        psd = df['psd'].tolist()
        ax.plot(time, psd, '-o', markersize= 4, color = colors[i], label = f'$K$ = {list_Ks[i]}')

    title = '$\\mu$ = ' + str(mu)
    ax.set_yscale('log')
    ax.set_ylabel(f'PSD [{unit_psd}]', fontsize=17)
    ax.set_xlabel('Time', fontsize=17)
    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
    ax.tick_params(axis='both', labelsize=16)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval = 1))  # Localizador de días
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))  # Formato de fecha
    ax.set_title(title, fontsize = 20)
    plt.show()


def psd_time_mus(list_dfs, list_mus, K, colors, unit_psd):

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)
    for i, df in enumerate(list_dfs):
        time = df['time'].tolist()
        psd = df['psd'].tolist()
        ax.plot(time, psd, '-o', markersize= 4, color = colors[i], label = f'$\\mu$ = {list_mus[i]}')

    title = f'$K$ = {str(K)}'
    ax.set_yscale('log')
    ax.set_ylabel(f'PSD [{unit_psd}]', fontsize=17)
    ax.set_xlabel('Time', fontsize=17)
    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
    ax.tick_params(axis='both', labelsize=16)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval = 1))  # Localizador de días
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))  # Formato de fecha
    ax.set_title(title, fontsize = 20)
    plt.show()
