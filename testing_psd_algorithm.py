#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Script para probar cómo calcular PSD
'''
from Download_data import download_ect as dd_ect
from Download_data import download_omni as dd_omni
from Process_data import process_ect as pp_ect
from Process_data import process_omni as pp_omni

import Source.inputs_MF_model as mf
import Source.plots_fedu_processing as plots_fp
import Source.fedu_processing as fp
import Source.invariants_calculation as inv
import Source.plots_invariants_calculation as plots_inv
import Source.plots_irbem_magnetic_field as plots_irb
import plots_poster_PSD_algorithm as poster

import os
import pathlib
import pandas as pd
from astropy import units as u
from astropy import constants as const
import matplotlib.pyplot as plt

import numpy as np
#import IRBEM

'''
################################################################################
###                 Definiciones preliminares del usuario
################################################################################
'''
############### Cosas que setea el usuario

### Modelo campo magnético a usar (opciones: T89, T96)
model = 'T89'

### Resolución datos omni a usar
resolution = "1h"
file_type = ''

### Promedio temporal para datos ECT
time_avg = '1min'

### Fechas análisis (YYYY-MM-DD)
sdate = '2013-05-30'
edate = '2013-05-30'

### Opciones para el modelo de campo mangético
# Queremos que calcule L*, el sistema de coordenadas input es 1:GEO
options = [1,0,0,0,0]

# definimos satélite a analizar ('a' o 'b')
sat = 'a'

### Valores de testeo para K, alphaK, mu y Emu.
# Fijar a priori un valor para alphaK y Emu no está bien (hay que modificarlo cuando sección 1 está lista)
target_mu = 0.02
target_K = 100

# Rangos de energía para los plots
E_min = 0.02*u.MeV # MeV
E_max = 15*u.MeV # MeV
energy_range_plot = [E_min, E_max]


### Cut offs canales de energía (para que no los usen)
# Rept tiene 12 canales de energía
rept_cut_offs = [1, 11]
# mageis tiene 25 canales de energía
mageis_cut_offs = [1, 14]

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

# Mostrar y guardar gráficos
show_flag = 1
save_flag = 0

'''
################################################################################
###                 Creamos constantes para el cálculo de PSD
################################################################################
'''

### Definimos constantes y unidades importantes
# Unidades de Re*(nT**(1/2)) para K
# Unidades de MeV/nT para mu
R_earth = 6371
unit_Re = u.def_unit('Re', const.R_earth.value*const.R_earth.unit)

unit_K = unit_Re*np.sqrt(1*u.nT)
target_K = target_K*unit_K

unit_mu = u.MeV/u.nT
target_mu = target_mu*unit_mu

unit_flux_mev =(u.cm**(2)*u.s*u.sr*u.MeV)**(-1)
unit_flux_kev =(u.cm**(2)*u.s*u.sr*u.keV)**(-1)

unit_c = u.def_unit('c', const.c)


'''
################################################################################
###       Definimos variables importantes para los inputs del MF model
################################################################################
'''

############### Definiciones importantes que no elige el usuario

### Fechas descarga omni
sdate_omni_1h, sdate_omni_min, edate_omni = mf.get_omni_dates(sdate, edate)

### Variables omni según resolución y modelo de campo
ect_var = mf.model_variables_ect(model)
omni_var = mf.model_variables_omni(model, resolution)
[var_ect, rename_ect] = ect_var
[var_omni_1h, var_omni_min, rename_omni] = omni_var

### Creamos el objeto de campo magnético T89
#model_obj = IRBEM.MagFields(options=options, kext=model, verbose = False, sysaxes=1)


'''
################################################################################
###                                 PARTE A:                                 ###
###     Descarga desde repositorio remoto: ECT-REP, ECT-MagEIS & OMNI        ###                                  ###
################################################################################
'''

############### Definiendo directorios importantes
cwd_path = pathlib.Path.cwd()

remote_dir_ect = "https://rbsp-ect.newmexicoconsortium.org/data_pub/"
remote_dir_omni = "https://cdaweb.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/"

data_dir_ect = str(cwd_path.parents[1]/'Data'/'rbsp'/'')
data_dir_omni = str(cwd_path.parents[1]/'Data'/'omni'/'')
plots_dir = str(cwd_path.parents[1]/'Plots'/'PSD_calculation'/'')

############### Descarga de datos ECT-REPT, ECT-MagEIS, y OMNI
'''
### REPT
dd_ect.download_CDFfiles_ECT(sdate, edate, remote_dir_ect, data_dir_ect,
                             probe=sat, instrument = 'rept', level="3")

### MagEIS
dd_ect.download_CDFfiles_ECT(sdate, edate, remote_dir_ect, data_dir_ect,
                             probe=sat, instrument = 'mageis', level="3")

### OMNI
dd_omni.download_CDFfiles_OMNI(sdate_omni_1h, edate_omni, remote_dir_omni,
                               data_dir_omni, '1h', '')

if resolution !='1h':
    dd_omni.download_CDFfiles_OMNI(sdate_omni_min, edate_omni, remote_dir_omni,
                                   data_dir_omni, resolution, file_type)
'''

'''
################################################################################
###                                 PARTE B:                                 ###
###     Carga desde repositorio local: ECT-REP, ECT-MagEIS & OMNI            ###                                 ###
################################################################################
'''
############### Lectura de datos ECT-REPT, ECT-MagEIS, y OMNI #################

### REPT
[rept_probe] = pp_ect.load_CDFfiles_ECT(sdate, edate, data_dir_ect,
               var_ect, rename_ect, sat, 'rept', '3', 'fedu')

[rept_info, rept_fedu, rept_info_meta, rept_fedu_meta] = rept_probe

### MagEIS
[mageis_probe] = pp_ect.load_CDFfiles_ECT(sdate, edate, data_dir_ect, var_ect,
                 rename_ect, sat, 'mageis', '3', 'fedu')

[mageis_info, mageis_fedu, mageis_info_meta, mageis_fedu_meta] = mageis_probe

### OMNI
omni_info_1h, omni_meta_1h = pp_omni.load_CDFfiles_OMNI(sdate_omni_1h,
                             edate_omni, data_dir_omni, var_omni_1h,
                             rename_omni, '1h', '')

omni_info_min = False
omni_meta_min = False

if resolution !='1h':
    omni_info_min, omni_meta_min = pp_omni.load_CDFfiles_OMNI(sdate_omni_min,
                                   edate_omni, data_dir_omni, var_omni_min,
                                   rename_omni, resolution, file_type)


'''
################################################################################
###                     PARTE C: Procesamiento de datos                      ###
###           Creación inputs para el modelo de campo magnético              ###
################################################################################
'''
### Juntamos info omni resolución de 1 hour & min
omni_info, omni_meta = mf.sync_omni_data([omni_info_1h, omni_meta_1h], [omni_info_min, omni_meta_min], resolution)

### Promediamos y sincronizamos
mageis_info, mageis_fedus = mf.time_average_ect_data(mageis_info, mageis_fedu, sdate, edate, time_avg)
rept_info, rept_fedus = mf.time_average_ect_data(rept_info, rept_fedu, sdate, edate, time_avg)

ect_info = mf.sync_ect_data(rept_info, mageis_info)

### Creamos los inputs
x_inputs, mag_inputs, N = mf.create_inputs_MF(ect_info, omni_info, R_earth, sdate, edate)

### Chequeamos que estén bien creados
mf.check_x_inputs(x_inputs, ect_info, R_earth)
mf.check_mag_inputs(mag_inputs, omni_info, ect_info, resolution)


'''
################################################################################
###                     PARTE D: Procesamiento datos                         ###
###           Agregamos unidades y fijamos un tiempo para los inputs         ###
################################################################################
'''
### Agregamos unidades de medida a energías, pitch angle y flujos
# Unidades para energy and alpha bins for REPT and MagEIS
rept_energy_bins =  rept_fedu_meta['energy_values']*u.MeV
rept_alpha_bins = rept_fedu_meta['alpha_values']*u.degree
rept_bins = [rept_energy_bins, rept_alpha_bins]

mageis_energy_bins =  mageis_fedu_meta['energy_values']*u.keV
mageis_alpha_bins = mageis_fedu_meta['alpha_values']*u.degree
mageis_bins = [mageis_energy_bins, mageis_alpha_bins]

# Unidades para fedus de REPT y MagEIS
rept_fedus_input = rept_fedus*unit_flux_mev
mageis_fedus_inputs = mageis_fedus*unit_flux_kev

# Elegimos un set de inputs para un tiempo determinado
#i = 180
i=0
x_input = x_inputs[i]
mag_input = mag_inputs[i]
inputs = [x_input, mag_input]
print(x_input)

rept_fedu = rept_fedus_input[i]
mageis_fedu = mageis_fedus_inputs[i]



'''
################################################################################
###                        PARTE E: Campo magnético                          ###
###                    Visualizando las líneas de campo                      ###
################################################################################
'''
### Graficamos la línea de campo que pasa por la posición del satélite
'''
R0 = 1
dens = 2
out_dict_fl = model_obj.trace_field_line(x_input, mag_input, R0)
posit_fl = out_dict_fl['POSIT']

ax = plots_irb.field_line(posit_fl, dens)
ax.scatter(x_input['x1'], x_input['x2'], x_input['x3'], color='red', s=100, label='initial pos')
ax.legend()

### Graficamos un drift shell de una partícula que está en la posición del satélite
# con pitch angle alpha
alpha = 60
out_dict_dbo = model_obj.drift_bounce_orbit(x_input, mag_input, alpha, R0)
posit = out_dict_dbo['POSIT']
Nposit = out_dict_dbo['Nposit']

dens = 3
title = '$\\alpha$ = 60°'
ax = plots_irb.drift_shell(posit, Nposit, dens, 2, title)
ax.scatter(x_input['x1'], x_input['x2'], x_input['x3'], color='r', s=100, label="Particle's position on " + str(x_input['dateTime'])[:-15] + 'at ' + str(x_input['dateTime'])[12:-7])
ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
plt.show()
'''

'''
################################################################################
###                         PARTE F: PSD calculation                         ###
###           Aplicamos algoritmo para calcular PSD for each timestamp       ###
################################################################################
'''


### Obtenemos info para los fits de pitch angle y energía
PA_fit_info = fp.info_fit_PA_flux(PA_fit_opt)
energy_fit_info = fp.info_fit_energy_data(energy_fit_opt)

# Obtenemos función para el cálculo de K
func_K, _, _ = inv.info_calculate_K(K_opt)


'''                    Section 1: Calculo de invariantes                 '''
'''
######## Step 1: Calcular K para distintos PA alphas, usando drift_bounce_orbit

# Array de PAs en grados
#alphas_lin = np.linspace(55, 90, 20)*u.degree
#alphas_geom = np.geomspace(3, 55, 20)*u.degree
alphas = np.linspace(3, 90, 40)*u.degree
#alphas_plot = np.linspace(0, 90, 100)*u.degree

# Obtenemos los valores de K para cada arreglo de PA en unidades de Re*(nT**(1/2))
#Ks_lin, _ = func_K(alphas_lin, model_obj, inputs, unit_K)
#Ks_geom, _ = func_K(alphas_geom, model_obj, inputs, unit_K)
Ks, _ = func_K(alphas, model_obj, inputs, unit_K)
#Ks_plot, _ = func_K(alphas_plot, model_obj, inputs, unit_K)


######## Step 2: Interpolamos la función alpha(K) y calculamos target alpha_K
# Limpiamos datos nan
#mask_lin = np.isnan(Ks_lin.value)
#Ks_lin = Ks_lin[~mask_lin]
#alphas_lin = alphas_lin[~mask_lin]

#mask_geom = np.isnan(Ks_geom.value)
#Ks_geom = Ks_geom[~mask_geom]
#alphas_geom = alphas_geom[~mask_geom]

mask = np.isnan(Ks.value)
Ks = Ks[~mask]
alphas = alphas[~mask]

# Obtenemos el interpolador spline para alpha(K)
#spline_lin = inv.interpolator_alpha(alphas_lin, Ks_lin)
#spline_geom = inv.interpolator_alpha(alphas_geom, Ks_geom)
spline  = inv.interpolator_alpha(alphas, Ks)

# Entre estos valores puedo elegir K
K_min = Ks[-1]
K_max = Ks[0]
#K_lim = Ks_lin[0]

#if target_K > K_lim:
#    spline_to_use = spline_geom
#else:
#    spline_to_use = spline_lin

# Obtenemos el valor de alpha para el target K
#target_alphaK = inv.interpolate_alpha_K(target_K, spline_to_use)
target_alphaK = inv.interpolate_alpha_K(target_K, spline)

# Visualizaciones
#list1_Ks = [Ks_plot, Ks_lin, Ks_geom]
#list2_Ks = [Ks_plot, Ks]

#list1_alphas = [alphas_plot, alphas_lin, alphas_geom]
#list2_alphas = [alphas_plot, alphas]

#plots_inv.interpolation_alpha_K(list1_Ks, list1_alphas, [spline_lin, spline_geom], K_min, K_max, K_lim)
#plots_inv.interpolation_alpha_K(list2_Ks, list2_alphas, [spline], K_min, K_max)
#poster.step2(list1_Ks, list1_alphas, [spline_lin, spline_geom], target_K, target_alphaK, K_min, K_max, K_lim, save=0)
#poster.step2(list2_Ks, list2_alphas, [spline], target_K, target_alphaK2, K_min, K_max)

#check_step2_180(Ks_lin, Ks_geom, model_t89, x_input, mag_input, Re, spline_lin, spline_geom, K_range)


######## Step 3: Calcular energía(mu,K) of chosen mu and K

# Obtenemos el campo magnético desde el modelo en la posición del spacecraft
# Esto lo tenemos que modificar para que use emfisis
dict_mag_field = model_obj.get_field_multi(x_input, mag_input)
b_mag = dict_mag_field['Bl'][0]*u.nT # Recordar que está en nanoteslas

# Entre estos valores puedo elegir mu
mu_min = inv.calculate_mu(E_min, b_mag, target_alphaK)
mu_max = inv.calculate_mu(E_max, b_mag, target_alphaK)

list_E = np.linspace(E_min, E_max, 100)
list_mu = inv.calculate_mu(list_E, b_mag, target_alphaK)

target_Emu = inv.calculate_E(target_mu, b_mag, target_alphaK)
#target_Emu2 = inv.calculate_E(target_mu, b_mag, target_alphaK2)

#plots_inv.calculation_E_mu(list_E, list_mu, target_Emu, target_mu)

######## Step 4: Calculate L*
func_Lstar, _, _ = inv.info_calculate_Lstar(Lstar_opt)
lstar = func_Lstar(target_alphaK, model_obj, inputs)

'''
'''                      Section 2: FEDU processing                      '''
target_alphaK = 38.48988795*u.degree
target_Emu = 2.75217844*u.MeV
PA_fit_opt = '2'
PA_fit_info = fp.info_fit_PA_flux(PA_fit_opt)
######## Step 5: Obtain j(alphaK) para cada canal de energía           #####

# Plot PA flux dist, j(alpha), para cada canal de energía de REPT y MagEIS
#plots_fp.PA_flux_data(rept_energy_bins, rept_alpha_bins, rept_fedu,
#                        'REPT', plots_dir, '', show_flag, save_flag)

#plots_fp.PA_flux_data(mageis_energy_bins, mageis_alpha_bins, mageis_fedu,
#                        'MagEIS', plots_dir, '', show_flag , save_flag)


### Ploteamos la función a ajustar, solo para visualizarla
#_, _ = plots_fp.PA_flux_function(*PA_fit_info[:-1], 1*unit_flux_mev, plots_dir,
#                                 show_flag, save_flag)


### Ajustar REPT y MagEIS PA flux y calculamos flux at target alphaK

rept_N = 12
mageis_N = 25

rept_channels = fp.channels_to_use(rept_cut_offs, rept_N)
mageis_channels = fp.channels_to_use(mageis_cut_offs, mageis_N)

''' Solo una función fit PA '''
print('ONE FUNCTION')
# Obtenemos el mejor ajuste, usando solo una función determinada de antes
rept_PA_fit_res = fp.fitting_alpha_flux(*PA_fit_info[2:], rept_fedu,
                  rept_channels, rept_alpha_bins, rept_N)
mageis_PA_fit_res = fp.fitting_alpha_flux(*PA_fit_info[2:], mageis_fedu,
                    mageis_channels, mageis_alpha_bins, mageis_N)

rept_fit, _, _, rept_err, rept_err_fit = rept_PA_fit_res
mageis_fit, _, _, mageis_err, mageis_err_fit = mageis_PA_fit_res

rept_fit_opts = [PA_fit_opt if val else None for val in rept_fit]
mageis_fit_opts = [PA_fit_opt if val else None for val in mageis_fit]

rept_func = [PA_fit_info[2] if val else None for val in rept_fit]
mageis_func = [PA_fit_info[2] if val else None for val in mageis_fit]

# Calculamos flujo at target alphaK
rept_flux_alphaK = fp.fitted_flux_at_alphaK(rept_PA_fit_res, rept_func, target_alphaK)
mageis_flux_alphaK = fp.fitted_flux_at_alphaK(mageis_PA_fit_res, mageis_func, target_alphaK)

flux_alphaK = [rept_flux_alphaK, mageis_flux_alphaK]
energy_bins = [rept_energy_bins[rept_fit], mageis_energy_bins[mageis_fit]]

### Checkeamos el fit

plots_fp.check_fit_PA_flux(rept_fit_opts, rept_func, rept_fedu, rept_bins, rept_PA_fit_res,
         rept_flux_alphaK, target_alphaK, 'REPT', plots_dir, show_flag, save_flag)


plots_fp.check_fit_PA_flux(mageis_fit_opts, mageis_func, mageis_fedu, mageis_bins,
         mageis_PA_fit_res, mageis_flux_alphaK, target_alphaK, 'MagEIS',
         plots_dir, show_flag, save_flag)

''' Selecciona la mejor función para fir PA '''
# Obtenemos el mejor ajuste, probando las 2 posibles funciones

print('BEST FUNCTION')
list_PA_fit_info = [fp.info_fit_PA_flux('1'), fp.info_fit_PA_flux('2')]

rept_PA_fit_res2 = fp.fitting_alpha_flux_V2(list_PA_fit_info, rept_fedu,
                   rept_channels, rept_alpha_bins, rept_N)

mageis_PA_fit_res2 = fp.fitting_alpha_flux_V2(list_PA_fit_info, mageis_fedu,
                     mageis_channels, mageis_alpha_bins, mageis_N)

rept_fit2, _, _, rept_err2, rept_err_fit2, rept_fit_opts2, rept_func2 = rept_PA_fit_res2
mageis_fit2, _, _, mageis_err2, mageis_err_fit2, mageis_fit_opts2, mageis_func2 = mageis_PA_fit_res2


rept_flux_alphaK2 = fp.fitted_flux_at_alphaK(rept_PA_fit_res2, rept_func2, target_alphaK)
mageis_flux_alphaK2 = fp.fitted_flux_at_alphaK(mageis_PA_fit_res2, mageis_func2, target_alphaK)

### Checkeamos el fit
plots_fp.check_fit_PA_flux(rept_fit_opts2, rept_func2, rept_fedu, rept_bins, rept_PA_fit_res2,
         rept_flux_alphaK2, target_alphaK, 'REPT', plots_dir, show_flag, save_flag)


plots_fp.check_fit_PA_flux(mageis_fit_opts2, mageis_func2, mageis_fedu, mageis_bins,
         mageis_PA_fit_res2, mageis_flux_alphaK2, target_alphaK, 'MagEIS',
         plots_dir, show_flag, save_flag)

#poster.step5(PA_fit_info, rept_fedu, rept_bins, rept_PA_fit_res, rept_flux_alphaK, target_alphaK, 0)
#poster.step5(PA_fit_info, mageis_fedu, mageis_bins, mageis_PA_fit_res, mageis_flux_alphaK, target_alphaK, 13)


'''

############### Step 6: Obtain psd(Emu) using j(alphaK) for each energy

### Mezclamos mageis y rept fluxes at target alphaK
energy_to_fit, energy_flux_to_fit = fp.join_energy_flux(flux_alphaK, energy_bins,
                                    unit_flux_mev)

### Calculamos psd to fit at target alphaK
energy_psd_to_fit = fp.flux_to_psd(energy_to_fit, energy_flux_to_fit[:, np.newaxis],
                                    unit_c, unit_flux_mev)
energy_psd_to_fit = energy_psd_to_fit[:,0]

### Plot energy flux j(E, alphaK) for target alphaK

flag_log = True
plots_fp.energy_y_data('flux', energy_to_fit, energy_flux_to_fit,
      target_alphaK, plots_dir, flag_log, show_flag, save_flag)

### Plot energy psd f(E, alphaK) for target alphaK
plots_fp.energy_y_data('psd', energy_to_fit, energy_psd_to_fit,
      target_alphaK, plots_dir, flag_log, show_flag, save_flag)


### Ajustar y=flux/psd at target alphaK y calculamos psd at target E_mu
Emu = np.array([target_Emu.value])*target_Emu.unit

if energy_fit_opt[1] == 'flux':
    y_data_to_fit = energy_flux_to_fit
    energy_fit_results = fp.fitted_y_at_Emu(energy_to_fit, y_data_to_fit,
                         *energy_fit_info[2:], 'flux', Emu)
    fit_obj, parms, flux_Emu = energy_fit_results
    psd_Emu = fp.flux_to_psd(Emu, flux_Emu, unit_c, unit_flux_mev)
elif energy_fit_opt[1] == 'psd':
    y_data_to_fit = energy_psd_to_fit
    energy_fit_results = fp.fitted_psd_at_Emu(energy_to_fit, y_data_to_fit,
                         *energy_fit_info[2:], 'psd', Emu)
    fit_obj, parms, psd_Emu = energy_fit_results
    flux_Emu = 0


plots_fp.check_fit_energy_data(energy_fit_results, energy_to_fit,
         y_data_to_fit, energy_fit_info, energy_range_plot, Emu,
         plots_dir, show_flag, save_flag)

#poster.step6(energy_fit_results, energy_to_fit, y_data_to_fit, energy_fit_info, energy_range_plot, Emu, save=0)
'''
