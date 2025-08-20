#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Script para calcular PSD para una serie de tiempo
'''

import os
import pathlib
import numpy as np
import Source.get_units as get_units
import Source.psd_processing as psd_proc
import Source.plots_psd as plots_psd
from itertools import product
import pandas as pd
'''
################################################################################
###                 Definiciones preliminares del usuario
################################################################################
'''

### Definiendo directorios importantes
cwd_path = pathlib.Path.cwd()
data_dir = str(cwd_path.parents[1]/'Data'/'PSD'/'')

### Modelo campo magnético a usar (opciones: T89, T96)
flag_model = 'T89'

### Promedio temporal para datos ECT
time_avg = '1min'

### Fechas análisis (YYYY-MM-DD HH:MM:SS)


# Definimos satélite a analizar ('a' o 'b')
sat = 'a'

### Valores para K y mu
# Unidades de Re*(nT**(1/2)) para K
# Unidades de MeV/nT para mu
units_inv = get_units.adiab_inv()
units_flux = get_units.flux()
units_psd = get_units.psd()

_, _, unit_K, unit_mu = units_inv
target_mus = [0.02, 0.01]*unit_mu
target_Ks = [60, 100, 200, 40]*unit_K
target_Ks = [100, 80]*unit_K

colors = ['b', 'g', 'r', 'k']

'''
################################################################################
###                              PARTE B: PSD                                ###
###           Aplicamos algoritmo para calcular PSD for each timestamp       ###
################################################################################
'''

'''
start_time = '2013-05-30 03:00:00'
end_time = '2013-05-30 07:00:00'

dict_dfs_Ks = {}
#for K, mu in product(target_Ks, target_mus):
for i, K in enumerate(target_Ks):
    print(i, K)
    list_mus = []
    for j, mu in enumerate(target_mus):
        print(j, mu)
        df, unit_psd = psd_proc.read_psd(data_dir, K, mu, start_time,
                                         flag_model, time_avg, 'storm')
        df_clean = psd_proc.clean_psd(df)
        list_mus.append(df_clean)
    dict_dfs_Ks[K] = list_mus

dict_dfs_mus = {}
#for K, mu in product(target_Ks, target_mus):
for i, mu in enumerate(target_mus):
    print(i, mu)
    list_Ks = []
    for j, K in enumerate(target_Ks):
        print(j, K)
        df, unit_psd = psd_proc.read_psd(data_dir, K, mu, start_time,
                                         flag_model, time_avg, 'storm')
        df_clean = psd_proc.clean_psd(df)
        list_Ks.append(df_clean)
    dict_dfs_mus[mu] = list_Ks


#plots_psd.psd_lstar(df, K, mu, unit_psd)
#plots_psd.psd_lstar(df_clean, K, mu, unit_psd)

#plots_psd.psd_time(df, K, mu, unit_psd)
#plots_psd.psd_time(df_clean, K, mu, unit_psd)


plots_psd.psd_lstar_Ks(dict_dfs_mus[mu], target_Ks, mu, colors, unit_psd)
plots_psd.psd_lstar_mus(dict_dfs_Ks[K], target_mus, K, colors, unit_psd)

plots_psd.psd_time_Ks(dict_dfs_mus[mu], target_Ks, mu, colors, unit_psd)
plots_psd.psd_time_mus(dict_dfs_Ks[K], target_mus, K, colors, unit_psd)
'''

'''
################################################################################
###                              PARTE B: PSD                                ###
###           Aplicamos algoritmo para calcular PSD for each timestamp       ###
################################################################################
'''
day = '2013-05-30'
times = ['03:00:00', '07:45:00', '12:15:00', '16:45:00']
delta = pd.Timedelta('4h')
K = target_Ks[0]
mu = target_mus[0]
df, unit_psd = psd_proc.read_psd(data_dir, K, mu, day, flag_model, time_avg, 'day')
df_clean = psd_proc.clean_psd(df)

dfs = []
for time in times:
    stime = pd.to_datetime(f'{day} {time}')
    etime = stime + delta

    mask = (df_clean['time'] >= stime) & (df_clean['time'] <= etime)
    df_filtered = df_clean[mask]
    dfs.append(df_filtered)
plots_psd.psd_lstar_times(dfs, times, mu, K, colors, unit_psd)
