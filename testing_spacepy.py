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
import spacepy
import spacepy.time
import spacepy.irbempy
import spacepy.coordinates

'''
################################################################################
###                 Definiciones preliminares del usuario
################################################################################
'''
### Flag para descargar datos
flag_dwnl= False

### Modelo campo magnético a usar (opciones: T89, T96)
flag_model = 'T89'

### Resolución datos omni a usar
res_omni = '1h'
file_type_omni = ''

### Promedio temporal para datos ECT
time_avg = '1min'

### Fechas análisis (YYYY-MM-DD HH:MM:SS)
start_time = '2013-05-30 00:00:00'
end_time = '2013-05-31 23:59:59'

start_time = '2013-05-30 09:00:00'
end_time = '2013-05-30 13:00:00'

sdate = start_time[:10]
edate = end_time[:10]

### Opciones para el modelo de campo mangético
# Queremos que calcule L*, el sistema de coordenadas input es 1:GEO
options = [1,0,0,0,0]

# Definimos satélite a analizar ('a' o 'b')
sat = 'a'

### Valores de testeo para K, alphaK, mu y Emu.
target_mus = [0.01, 0.02]
target_Ks = [100, 80]

# Rangos de energía para los plots
E_min = 0.01*u.MeV # MeV
E_max = 15*u.MeV # MeV
energy_range_plot = [E_min, E_max]


### Cut offs canales de energía (para que no los usen)
# Rept tiene 12 canales de energía
rept_N = 12
rept_cut_offs = [1, 12]
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



'''
################################################################################
###                      PARTE B: INPUTS IRBEM Y SPACEPY                                ###
###     Creación inputs para los modelos de campo magnético y cálculo psd      ###
################################################################################
'''

### Inputs para irbem
inputs_mf, inputs_fedus, N = psd.step0_inputs(ect_info, omni_info, start_time,
                             end_time, res_omni, omni_range, rept_fedus, mageis_fedus,
                             units)
x_inputs, mag_inputs = inputs_mf

### Inputs para spacepy
# Datos ect
ect_filt = mf.filter_dates(ect_info, start_time, end_time)
idx1 = ect_filt.index[0]
idx2 = ect_filt.index[-1]
dates, xs_geo, ys_geo, zs_geo = mf.get_GEO_coordinates(ect_filt, units_inv[0])

times = spacepy.time.Ticktock(list(dates))
list_coord = [[a, b, c] for a, b, c in zip(xs_geo, ys_geo, zs_geo)]
coords = spacepy.coordinates.Coords(list_coord, dtype='GEO', carsph='car')

'''
'''
################################################################################
###                   PARTE C: Cálculo invariantes usando spacepy            ###
###           Usamos la lista completa, no recorremos cada elemento          ###
################################################################################

magkeys = {'Kp':0, 'Dst':1, 'dens':2, 'velo':3, 'Pdyn':4, 'ByIMF':5, 'BzIMF':6,
           'G1':7, 'G2':8, 'G3':9, 'W1':10, 'W2':11, 'W3':12, 'W4':13, 'W5':14,
           'W6':15}

# Array de PAs en grados
alphas_right = np.geomspace(30, 90, 21)[1:]*u.degree
alphas_left = np.geomspace(30, 3, 20)[::-1]*u.degree
alphas_total = np.concatenate([alphas_left, alphas_right])
#spacepy.config['ncpus'] = 1
n = spacepy.config['ncpus']
print(n)
# Step 1: calculamos K para el arreglo de alphas

if n ==1:
    dictK_left, args_left = spacepy.irbempy.get_Lstar(times, coords, extMag = flag_model, options = options, alpha = alphas_left)
    dictK_right, args_right = spacepy.irbempy.get_Lstar(times, coords, extMag = flag_model, options = options, alpha = alphas_right)

else:
    dictK_left = spacepy.irbempy.get_Lstar(times, coords, extMag = flag_model, options = options, alpha = alphas_left)
    dictK_right = spacepy.irbempy.get_Lstar(times, coords, extMag = flag_model, options = options, alpha = alphas_right)

bmirr_left = dictK_left['Bmirr']
I_left = dictK_left['Xj']
K_left = I_left*np.sqrt(bmirr_left)

bmirr_right = dictK_right['Bmirr']
I_right = dictK_right['Xj']
K_right = I_right*np.sqrt(bmirr_right)

Ks = np.concatenate((K_left, K_right), axis=1)*unit_K

### calculamos el campo magnético
dictB = spacepy.irbempy.get_Bfield(times, coords, extMag=flag_model, options=options)
blocal = dictB['Blocal']*u.nT

if n ==1:
    magin_spy_left = {}
    magin_spy_right = {}
    for key in omni_range:
        magin_spy_left[key] = args_left[magkeys[key]]
        magin_spy_right[key] = args_right[magkeys[key]]

alphasK = []
Esmu = []
Lstar = []
for i in range(len(times)):
    list_K = Ks[i]
    mask = np.isnan(list_K)
    K = list_K[~mask]
    alphas = alphas_total[~mask]
    target_alphasK, _, _, _ = psd.step2(alphas, K, target_Ks)
    alphasK.append(target_alphasK)
    b = blocal[i]
    target_Esmu = psd.step3(b, target_alphasK, target_mus)
    Esmu.append(target_Esmu)

    dictL, _ = spacepy.irbempy.get_Lstar(times[i], coords[i], extMag = flag_model,
            options = options, alpha = target_alphasK)
    Lstar.append(dictL['Lstar'])


#mf.check_inputs_spy(mag_inputs, magin_spy_left, magin_spy_right, len(times), list(omni_range.keys()))

'''
################################################################################
###                  PARTE D: Cálculo invariantes usando IRBEM               ###
###                    Iteramos por cada elemento de la lista                ###
################################################################################
'''

### Obtenemos función para el cálculo de K y Lstar
info_K = inv.info_calculate_K(K_opt)
info_Lstar = inv.info_calculate_Lstar(Lstar_opt)

# Creamos los dataframes para guardar datos:
labels_Ks = [str(x) for x in target_Ks]
labels_mus = [str(x) for x in target_mus]

cols = pd.MultiIndex.from_product([labels_Ks, labels_mus], names=['K', 'mu'])
df_psd = pd.DataFrame(columns = cols)
df_lstar = pd.DataFrame(columns=cols)


# Loop en el tiempo
N1 = 50
N2 = 51
show = True
for i in range(N1, N2):
    print(i)
    # Obtenemos los inputs del modelo de campo magnético, para el timestep i
    inputs = [x_inputs[i], mag_inputs[i]]
    timestamp = x_inputs[i]['dateTime']
    print(f'TIMESTAMP {i}: ', timestamp)
    save_lstar = [df_lstar, timestamp, target_Ks, target_mus]


    ###################  Section 1: Calculo de invariantes  ####################

    ''' Step 1: Calcular K para distintos PA alphas '''
    ### Obtenemos los valores de K para el arreglo de PA en unidades de Re*(nT**(1/2))
    alphas, Ks = psd.step1(model_obj, inputs, alphas_total, info_K, unit_K)


    ''' Step 2: Interpolar la función alpha(K) y calcular target alpha_K '''
    target_alphasK, spline, K_min, K_max = psd.step2(alphas, Ks, target_Ks)
    target_alphasK = [38, 58, ]*u.deg
    ### Visualización
#    if show:
#        plots_psd_calc.check_step2(info_K, model_obj, inputs, unit_K, Ks, alphas,
#                              spline, target_Ks, target_alphasK)

    ''' Step 3: Calcular energía(mu,K) of chosen mu and K '''
    b_mag = inv.get_B_model(model_obj, inputs)
    target_Esmu = psd.step3(b_mag, target_alphasK, target_mus)

    ''' Step 4: Calculate L* '''
    df_lstar = psd.step4(model_obj, inputs, target_alphasK, info_Lstar,
                         save_lstar)
