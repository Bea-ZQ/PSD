#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Script para calcular PSD para una serie de tiempo
'''

import Source.inputs_MF_model as mf
import Source.fedu_processing as fp
import Source.invariants_calculation as inv
import Source.psd_calculation as psd
import Source.psd_processing as psd_proc
import Source.get_units as get_units

import Source.plots_psd_calculation as plots_psd_calc

import os
import pathlib
import pandas as pd
from astropy import units as u
import matplotlib.pyplot as plt

import numpy as np
import IRBEM

'''
################################################################################
###                 Definiciones preliminares del usuario
################################################################################
'''
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
start_time = '2013-05-30 00:00:00'
end_time = '2013-05-31 23:59:59'

#start_time = '2013-05-30 03:00:00'
#end_time = '2013-05-30 07:00:00'

sdate = start_time[:10]
edate = end_time[:10]

### Opciones para el modelo de campo mangético
# Queremos que calcule L*, el sistema de coordenadas input es 1:GEO
options = [1,0,0,0,0]

# Definimos satélite a analizar ('a' o 'b')
sat = 'a'

### Valores de testeo para K, alphaK, mu y Emu.
target_mus = [0.02, 0.01]
target_Ks = [100, 80, 120]
#target_Ks = [60]

# Rangos de energía para los plots
E_min = 0.01*u.MeV # MeV
E_max = 10*u.MeV # MeV
energy_range_plot = [E_min, E_max]


### Cut offs canales de energía (para que no los usen)
# Rept tiene 12 canales de energía
rept_N = 12
rept_cut_offs = [1, 6]
rept_channels = fp.channels_to_use(rept_cut_offs, rept_N)

# mageis tiene 25 canales de energía
mageis_N = 25
mageis_cut_offs = [1, 14]
mageis_channels = fp.channels_to_use(mageis_cut_offs, mageis_N)


### Fijamos las opciones para los fits de los flujos
# Pitch angle (opciones: '1' y '2')
PA_fit_opt = '2'
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
###                 Creamos constantes para el cálculo de PSD
################################################################################
'''

### Definimos unidades importantes
# Unidades de Re*(nT**(1/2)) para K
# Unidades de MeV/nT para mu

units_inv = get_units.adiab_inv()
units_flux = get_units.flux()
units_psd = get_units.psd()

_, _, unit_K, unit_mu = units_inv
units = [units_inv, units_flux, units_psd]

unit_flux_mev, unit_flux_kev = units_flux

target_Ks = target_Ks*unit_K
target_mus = target_mus*unit_mu

### Creamos el objeto de campo magnético T89
model_obj = IRBEM.MagFields(options=options, kext=flag_model, verbose = False, sysaxes=1)


'''
################################################################################
###                             PARTE A: DATA                                ###
###     Descarga desde repositorio remoto: ECT-REPT, ECT-MagEIS & OMNI       ###
###     Carga desde repositorio local: ECT-REPT, ECT-MagEIS & OMNI           ###
###     Creación inputs para el modelo de campo magnético y cálculo psd      ###
################################################################################
'''

### Definiendo directorios importantes
cwd_path = pathlib.Path.cwd()

data_dir_ect = str(cwd_path.parents[1]/'Data'/'rbsp'/'')
data_dir_omni = str(cwd_path.parents[1]/'Data'/'omni'/'')

### Descarga y lectura de datos ECT-REPT, ECT-MagEIS, y OMNI
### ECT
ect_info, rept, mageis = psd.step0_get_ECT_data(start_time, end_time, sat,
                         data_dir_ect, flag_dwnl, flag_model, time_avg)

rept_fedus, rept_bins = rept
mageis_fedus, mageis_bins = mageis

### OMNI
omni_info, omni_meta, omni_range = psd.step0_get_OMNI_data(start_time, end_time, data_dir_omni,
                       res_omni, file_type_omni, flag_dwnl, flag_model)

### Creamos los inputs
inputs_mf, inputs_fedus, N = psd.step0_inputs(ect_info, omni_info, start_time,
                             end_time, res_omni, omni_range, rept_fedus, mageis_fedus,
                             units)

x_inputs, mag_inputs = inputs_mf
rept_fedus_inputs, mageis_fedus_inputs = inputs_fedus


'''
################################################################################
###                              PARTE B: PSD                                ###
###           Aplicamos algoritmo para calcular PSD for each timestamp       ###
################################################################################
'''

### Obtenemos info para los fits de pitch angle y energía. Es siempre la misma
PA_fit_info = fp.info_fit_PA_flux(PA_fit_opt)
energy_fit_info = fp.info_fit_energy_data(energy_fit_opt)

### Obtenemos función para el cálculo de K y Lstar
info_K = inv.info_calculate_K(K_opt)
info_Lstar = inv.info_calculate_Lstar(Lstar_opt)

# Fijamos cuál función para el step6 vamos a usar, dependiendo de si queremos
# interpolar f(E) o j(E)
if energy_fit_opt[1] == 'flux':
    step6 = psd.step6_flux
elif energy_fit_opt[1] == 'psd':
    step6 = psd.step6_psd

# Array de PAs en grados
alphas_right = np.geomspace(30, 90, 21)[1:]*u.degree
alphas_left = np.geomspace(30, 3, 20)[::-1]*u.degree
alphas_total = np.concatenate([alphas_left, alphas_right])

# Creamos los dataframes para guardar datos:
labels_Ks = [str(x) for x in target_Ks]
labels_mus = [str(x) for x in target_mus]

#cols_lstar = pd.Index(labels_Ks, name='K')

cols = pd.MultiIndex.from_product([labels_Ks, labels_mus], names=['K', 'mu'])
df_psd = pd.DataFrame(columns = cols)
df_lstar = pd.DataFrame(columns=cols)


# Creamos los dataframes para guardar errores

rept_cols = [f'REPT_{channel}' for channel in range(1, rept_N+1)]
mageis_cols = [f'MagEIS_{channel}' for channel in range(1, mageis_N+1)]
all_cols = rept_cols + mageis_cols

df_rmse = pd.DataFrame(columns=all_cols)
df_r2 = pd.DataFrame(columns=all_cols)

df_rmse2 = pd.DataFrame(columns=all_cols)
df_r22 = pd.DataFrame(columns=all_cols)

# Loop en el tiempo
N = 2
show = False
for i in range(N):
    # Obtenemos los inputs del modelo de campo magnético, para el timestep i
    inputs = [x_inputs[i], mag_inputs[i]]
    timestamp = x_inputs[i]['dateTime']
    print(f'TIMESTAMP {i}: ', timestamp)
    # Obtenemos los flujo de REPT y MagEIS para el timestep i
    rept_fedu = rept_fedus_inputs[i]
    mageis_fedu = mageis_fedus_inputs[i]

    save_psd = [df_psd, timestamp, target_Ks]
    save_lstar = [df_lstar, timestamp, target_Ks, target_mus]


    ###################  Section 1: Calculo de invariantes  ####################

    ''' Step 1: Calcular K para distintos PA alphas '''
    ### Obtenemos los valores de K para el arreglo de PA en unidades de Re*(nT**(1/2))
    alphas, Ks = psd.step1(model_obj, inputs, alphas_total, info_K, unit_K)

    ''' Step 2: Interpolar la función alpha(K) y calcular target alpha_K '''
    target_alphasK, spline, K_min, K_max = psd.step2(alphas, Ks, target_Ks)
#    target_alphasK = [38, 58, ]*u.deg
    ### Visualización
    if show:
        plots_psd_calc.check_step2(info_K, model_obj, inputs, unit_K, Ks, alphas,
                              spline, target_Ks, target_alphasK)

    ''' Step 3: Calcular energía(mu,K) of chosen mu and K '''
    target_Esmu = psd.step3_model(model_obj, inputs, target_alphasK, target_mus,
                  energy_range_plot)
#    target_Esmu = [[0.5]*u.MeV, [2]*u.MeV]

    ''' Step 4: Calculate L* '''
    df_lstar = psd.step4(model_obj, inputs, target_alphasK, info_Lstar,
                         save_lstar)


    #######################  Section 2: FEDU processing  #######################

    ''' Step 5: Obtain j(alphaK) for each energy channel '''

    ''' using only one function'''
    save_errors = [df_rmse, df_r2, timestamp]
    ### Ajustar REPT y MagEIS PA flux y calculamos flux at target alphaK
    rept_flux_alphasK, rept_energy, rept_PA_fit_res= psd.step5_onefunc(PA_fit_info, rept_fedu,rept_channels, rept_bins, rept_N, target_alphasK, 'REPT', save_errors)
    mageis_flux_alphasK, mageis_energy, mageis_PA_fit_res = psd.step5_onefunc(PA_fit_info, mageis_fedu, mageis_channels, mageis_bins, mageis_N, target_alphasK, 'MagEIS', save_errors)

    rept_fit_opts = [PA_fit_opt if val else None for val in rept_PA_fit_res[0]]
    mageis_fit_opts = [PA_fit_opt if val else None for val in mageis_PA_fit_res[0]]

    rept_func = [PA_fit_info[2] if val else None for val in rept_PA_fit_res[0]]
    mageis_func = [PA_fit_info[2] if val else None for val in mageis_PA_fit_res[0]]

    ### Checkeamos el fit
    if show:
        plots_psd_calc.check_step5(rept_fit_opts, rept_func, rept_fedu, rept_bins,
                  rept_PA_fit_res, rept_flux_alphasK, target_alphasK, 'REPT')

        plots_psd_calc.check_step5(mageis_fit_opts, mageis_func, mageis_fedu, mageis_bins,
                  mageis_PA_fit_res, mageis_flux_alphasK, target_alphasK,
                  'MagEIS')

    flux_alphasK = [rept_flux_alphasK, mageis_flux_alphasK]
    energy_bins = [rept_energy, mageis_energy]


    #rmse, r2
    save_errors2 = [df_rmse2, df_r22, timestamp]
    ''' Selecting best function'''
    list_PA_fit_info = [fp.info_fit_PA_flux('1'), fp.info_fit_PA_flux('2')]

    rept_flux_alphasK2, rept_energy2, rept_PA_fit_res2 = psd.step5_bestfunc(list_PA_fit_info, rept_fedu,rept_channels, rept_bins, rept_N, target_alphasK, 'REPT', save_errors2)
    mageis_flux_alphasK2, mageis_energy2, mageis_PA_fit_res2 = psd.step5_bestfunc(list_PA_fit_info, mageis_fedu, mageis_channels, mageis_bins, mageis_N, target_alphasK, 'MagEIS', save_errors2)

    ### Checkeamos el fit
    if show:
        plots_psd_calc.check_step5(rept_PA_fit_res2[5], rept_PA_fit_res2[6], rept_fedu, rept_bins,
                  rept_PA_fit_res2, rept_flux_alphasK2, target_alphasK, 'REPT')

        plots_psd_calc.check_step5(mageis_PA_fit_res2[5], mageis_PA_fit_res2[6], mageis_fedu, mageis_bins,
                  mageis_PA_fit_res2, mageis_flux_alphasK2, target_alphasK,
                  'MagEIS')


    ''' Step 6: Obtain psd(Emu) using j(alphaK) for each energy '''

    ### Mezclamos mageis y rept fluxes at target alphaK
    energy_to_fit, energy_flux_to_fit = fp.join_energy_flux(flux_alphasK,
                                        energy_bins, unit_flux_mev)

    ### Ajustar y=flux/psd at target alphaK y calculamos psd at target E_mu
    df_psd, fit_objs, pars = step6(energy_to_fit, energy_flux_to_fit,
                             energy_fit_info, target_Esmu, target_alphasK,
                             units_psd, save_psd)

    ### Checkeamos los datos y el fit
    if show:
        plots_psd_calc.check_step6(energy_flux_to_fit, energy_to_fit, target_alphasK,
        target_Ks, fit_objs, pars, df_psd, energy_fit_info, energy_range_plot,
        target_Esmu, units_psd, fp.flux_to_psd)
    show = False



plots_psd_calc.psd_lstar(df_psd, df_lstar, target_Ks[0], target_mus[0])
plots_psd_calc.psd_time(df_psd, target_Ks[0], target_mus[0])
