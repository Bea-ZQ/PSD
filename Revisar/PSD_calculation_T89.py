###############################################################################
### En este script estoy calculando invariantes adiabáticas.
################################################################################
import os, sys
import datetime
import pandas as pd
import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as u
import IRBEM

from PSD_calculation_utils import step0_get_RETP_data, step0_get_OMNI_data, step0_get_inputs
from PSD_calculation_utils import step1_get_Ki_alphai, step2_get_target_alphaK
from PSD_calculation_utils import step3_get_target_Emu, step4_get_j_alphaK
from PSD_calculation_utils import step5_6_get_f_Emu, step7_get_Lstar
from PSD_calculation_utils import checking_step2, checking_step3, checking_step4
from PSD_calculation_utils import checking_step5_6, plot_psd_lstar, plot_psd_time, plot_flux_time


################################################################################
###  Definiciones preliminares
################################################################################

### Unidades de Re*(nT**(1/2)) para K
### Unidades de MeV/nT para mu

### Definimos target K y mu para el cálculo
target_Ks = [100]*K_unit
#target_Ks=[415, 200]*K_unit
target_mu = 0.02*mu_unit

### Periodo de tiempo a analizar
start_time1 = sdate + ' ' + '03:00:00'
start_time2 = sdate + ' ' + '07:45:00'
start_time3 = sdate + ' ' + '12:15:00'
start_time4 = sdate + ' ' + '16:45:00'
start_time5 = edate + ' ' + '06:00:00'
start_time6 = edate + ' ' + '15:00:00'
d_hour = 4
start_times = [start_time1, start_time2, start_time3]
start_times = [start_time1, start_time2]
start_times = [start_time1, start_time5]
start_times = [start_time1, start_time2, start_time3, start_time4, start_time5, start_time6]
start_times = [start_time1, start_time2]

step = 15

################################################################################
###  Step 0: Get OMNI and REPT data between sdate and edate. Only return
###  relevant keys for magnetic field model
################################################################################

### Get omni data
parameters_T89 = ['Epoch', 'YR', 'Day', 'HR', 'KP']
rename_mapping = {'KP': 'Kp', 'Epoch': 'Epoch', 'YR': 'year', 'Day':'day', 'HR': 'hour'}
omni1h_info, omni1h_metadata = step0_get_OMNI_data(sdate, edate, False,parameters_T89, rename_mapping)

### Get REPT data
probeA_data = step0_get_RETP_data(sdate, edate, False, 'a')
#probeB_data = step0_get_RETP_data(sdate, edate, False, 'b')

[reptA_info, feduA, reptA_metadata, feduA_metadata] = probeA_data[0]
#[reptB_info, feduB, reptB_metadata, feduB_metadata] = probeB_data[0]

colors = ['r', 'b', 'g', 'cyan', 'yellow', 'k']

psd_data = []
flux_data = []
time_data = []
lstar_data = []
for j in range(len(start_times)):
    print('START TIME', start_times[j])

    '''
    ###  Step 0: Get x and magnetic inputs for analysis.
    '''

    list_x_input, list_mag_input, list_fedu_input, fedu_metadata, len = step0_get_inputs(start_times[j], d_hour, reptA_info, feduA, omni1h_info, feduA_metadata, R_earth)
    #list_x_input, list_mag_input, list_fedu_input, fedu_metadata, len = step0_get_inputs(start_time, d_hour, reptB_info, feduB, omni1h_info, feduB_metadata, R_earth)



    '''
    ### PSD transformation algorithm
    '''
#    len =1
    print(len)

    list_lstar = []
    list_psd = []
    list_flux = []
    list_time = []
    for i in range(0, len, step):
        print('Time step: ', i)
        x_input = list_x_input[i]
        mag_input = list_mag_input[i]
        model_info = [model_t89, x_input, mag_input]
        fedu_input = list_fedu_input[i]
        list_time.append(x_input['dateTime'])
        '''
        ### Step 1: Calculamos K para distintos alphas usando drift_bounce_orbit
        '''

        ### Calculamos Ki(alphai)
        data_lin, data_geom = step1_get_Ki_alphai(model_info, '2', Re_unit)
        [alphas_lin, Ks_lin] = data_lin
        [alphas_geom, Ks_geom] = data_geom

        '''
        ### Step 2: Interpolamos la función alpha(K) para obtener target alpha_K
        '''

        target_alphasK, spl_lin, spl_geom = step2_get_target_alphaK(data_lin, data_geom, target_Ks, flag180 = False)
#        checking_step2(Ks_lin, Ks_geom, model_info, '2', Re_unit, spl_lin, spl_geom, target_Ks, target_alphasK, False)
#        checking_step2(Ks_lin, Ks_geom, model_info, '2', Re_unit, spl_lin, spl_geom, target_Ks, target_alphasK, True)

        '''
        # Step 3: Calcular target Emu. Energía of chosen mu and K
        '''

        target_Esmu, b_mag  = step3_get_target_Emu('model', model_info, target_mu, target_alphasK)
    #    checking_step3(b_mag, target_alphasK, target_Esmu, target_mu)

        '''
        ### Step 4: Interpolar j(alpha) para cada canal de energía y calcular j(alphaK)
        '''

        ### Agregamos unidades a los datos

        #energy 1,8 MeV, alpha = 68°
        fedu_energy0_alpha6 = fedu[6,0].value
        list_flux.append(fedu_energy0_alpha6)
    #    js_alphaK = step4_get_j_alphaK(fedu, energy_bins, alpha_bins, '1', cut_off_index,
    #                                   target_alphasK, False)
        js_alphasK, bins_to_use, parameters, max_values  = step4_get_j_alphaK(fedu, energy_bins, alpha_bins, '2', cut_off_index, target_alphasK, False)
    #    checking_step4(fedu, alpha_bins, energy_bins,bins_to_use, parameters, max_values, '2', target_alphasK, js_alphasK)

        '''
        ### Step 5: Interpolar j(E) or f(E) using calculated fluxes/PSD at target alphaK y calcular j(Emu) or f(Emu)
        ### Step 6: Transformation from fluxes to PSD
        '''


        psd_Esmu_pl, par_opt_pl = step5_6_get_f_Emu(js_alphasK, energy_bins, bins_to_use, 'flux', 'pl', target_alphasK, target_Esmu, c_unit, False)
        #psd_Esmu_exp, par_opt_exp = step5_6_get_f_Emu(js_alphasK, energy_bins, bins_to_use, 'psd', 'exp', target_alphasK, target_Esmu, c_unit, False)
        #psd_Esmu_spl, par_opt_spl = step5_6_get_f_Emu(js_alphasK, energy_bins,bins_to_use, 'flux', 'spl', target_alphasK, target_Esmu, c_unit, False)

#        checking_step5_6(energy_bins, bins_to_use, 'flux', par_opt_pl, 'pl', js_alphasK, target_Esmu, psd_Esmu_pl, c_unit)
        #checking_step5_6(energy_bins, bins_to_use,'psd', par_opt_exp, 'exp', js_alphasK, target_Esmu, psd_Esmu_exp, c_unit)
        #checking_step5_6(energy_bins, bins_to_use,'flux', par_opt_spl, 'spl', js_alphasK, target_Esmu, psd_Esmu_spl, c_unit)
        list_psd.append(psd_Esmu_pl[0].value)


        '''
        ### Step 7: Calculate L*
        '''
        ls_star = step7_get_Lstar(model_info, '1', target_alphasK)
        list_lstar.append(ls_star[0])

    psd = list_psd*psd_Esmu_pl.unit
    flux = list_flux *fedu.unit

    psd_data.append(psd)
    flux_data.append(flux)
    time_data.append(list_time)
    lstar_data.append(list_lstar)

colors = ['r', 'b', 'darkturquoise', 'darkgreen', 'darkorange', 'k']
plot_psd_lstar(lstar_data, psd_data, start_times, colors, target_Ks[0], target_mu)
plot_psd_time(time_data, psd_data, start_times, colors, target_Ks[0], target_mu)
plot_flux_time(time_data, flux_data, start_times, colors, energy_bins[0], alpha_bins[6])
# Por ahora funciona solo con un K y un mu
