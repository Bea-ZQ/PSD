#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Script para calcular PSD para una serie de tiempo
'''

import os
import pathlib
import numpy as np
import Source.get_units as get_units
import Source.psd_calculation as psd_calc
import Source.psd_processing as psd_proc
import Source.plots_psd_calculation as plots_psd_calc

'''
################################################################################
###                 Definiciones preliminares del usuario
################################################################################
'''

### Definiendo directorios importantes
cwd_path = pathlib.Path.cwd()
data_dir_ect = str(cwd_path.parents[1]/'Data'/'rbsp'/'')
data_dir_omni = str(cwd_path.parents[1]/'Data'/'omni'/'')
save_dir = str(cwd_path.parents[1]/'Data'/'PSD'/'')
### Flag para descargar datos
flag_dwnl= False

### Modelo campo magnético a usar (opciones: T89, T96)
flag_model = 'T96'

### Resolución datos omni a usar
res_omni = '1h'
file_type_omni = ''

### Promedio temporal para datos ECT
time_avg = '1min'

### Fechas análisis (YYYY-MM-DD HH:MM:SS)
#start_time = '2013-02-28 18:15:00'
#end_time = '2013-03-03 11:23:00'

#start_time = '2013-05-30 00:00:00'
#end_time = '2013-05-30 23:59:59'
#type = 'day'

start_time = '2013-05-30 03:00:00'
end_time = '2013-05-30 07:00:00'
type = 'storm'

### Opciones para el modelo de campo (calcula L* y en coordenadas 1:GEO)
opt_mf = [1,0,0,0,0]

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

# Rangos de energía para los plots
E_min = 0.01 # MeV
E_max = 10 # MeV

### Cut offs canales de energía (para que no los usen)
# REPT tiene 12 canales de energía
rept_N = 12
rept_cut_offs = [1, 6]

# MagEIS tiene 25 canales de energía
mageis_N = 25
mageis_cut_offs = [1, 14]

### Fijamos las opciones para los fits de los flujos
# Pitch angle (opciones: '1' y '2')
PA_fit_opt = '2'
PA_fit_opt = ['1', '2']

# Energía (opciones: ('pl', 'flux'), ('exp', 'psd'), ('spl', 'flux')
#                    ('spl', 'psd'), ('lin', 'flux'), ('lin', 'psd'))
energy_fit_opt = ('lin', 'flux')

### Fijamos la opción para el cálculo de K
# Opciones (1: mirror points, 2: drift_bounce, 3:integración propia)
K_opt = '2'

### Fijamos la opción para el cálculo de Lstar
# Opciones (1: mirror points, 2: drift_bounce)
Lstar_opt = '2'


'''
################################################################################
###                             PARTE A: DATA                                ###
###     Descarga desde repositorio remoto: ECT-REPT, ECT-MagEIS & OMNI       ###
###     Carga desde repositorio local: ECT-REPT, ECT-MagEIS & OMNI           ###
###     Creación inputs para el modelo de campo magnético y cálculo psd      ###
################################################################################
'''
### Encapsulamos las variables para los cálculos psd
opts_psd = [PA_fit_opt, energy_fit_opt, K_opt, Lstar_opt]
opts_model = [opt_mf, flag_model]
units = [units_inv, units_flux, units_psd]
channels_to_use = [rept_cut_offs, rept_N, mageis_cut_offs, mageis_N]
targets = [target_Ks, target_mus]
times = [start_time, end_time]
data_dirs = [data_dir_ect, data_dir_omni]
energy_range = [E_min, E_max]

### Obtenemos los inputs
inputs, bins = psd_calc.get_inputs(times, sat, data_dirs, flag_dwnl, opts_model,
               time_avg, res_omni, file_type_omni, targets, units)


'''
################################################################################
###                              PARTE B: PSD                                ###
###           Aplicamos algoritmo para calcular PSD for each timestamp       ###
################################################################################
'''

df_psd, df_lstar, df_rmse, df_r2 = psd_calc.psd_calculation(channels_to_use, opts_psd, opts_model,
                   targets, inputs, energy_range, bins, units, N_steps = 10e9)



plots_psd_calc.psd_lstar(df_psd, df_lstar, target_Ks[0], target_mus[0])
plots_psd_calc.psd_time(df_psd, target_Ks[0], target_mus[0])


psd_proc.save_psd_data(save_dir, df_psd, df_lstar, df_rmse, df_r2, targets, times,
                  flag_model, time_avg, sat, type)
