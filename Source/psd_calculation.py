#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
import contextlib
from astropy import units as u
import pandas as pd
import IRBEM

from Download_data import download_ect as dd_ect
from Download_data import download_omni as dd_omni
from Process_data import process_ect as pp_ect
from Process_data import process_omni as pp_omni

import Source.invariants_calculation as inv
import Source.fedu_processing as fp
import Source.inputs_MF_model as mf


def step0_get_ECT_data(stime, etime, probe_opt, data_dir, flag_download,
                       flag_model, time_avg):
    print('------------------------------------------------------------------------')
    print(f'STEP 0: Get and process ECT (REPT & MagEIS) data for probe {probe_opt}')
    print(f'    Only return relevant keys for magnetic field model {flag_model}')
    print('------------------------------------------------------------------------')

    ### Fechas descarga ect
    sdate_ect = stime[:10]
    edate_ect = etime[:10]

    if flag_download:
        remote_dir = "https://rbsp-ect.newmexicoconsortium.org/data_pub/"

        # REPT
        dd_ect.download_CDFfiles_ECT(sdate_ect, edate_ect, remote_dir, data_dir,
            probe=probe_opt, instrument = 'rept', level="3")

        # MagEIS
        dd_ect.download_CDFfiles_ECT(sdate_ect, edate_ect, remote_dir, data_dir,
                probe=probe_opt, instrument = 'mageis', level="3")

    else:
        print('* No download')

    ### Variables ECT según modelo de campo magnético
    ect_var, ect_rename = mf.model_variables_ect(flag_model)

    ### REPT
    [rept_probe] = pp_ect.load_CDFfiles_ECT(sdate_ect, edate_ect, data_dir,
                   ect_var, ect_rename, probe_opt, 'rept', '3', 'fedu')

    ### MagEIS
    [mageis_probe] = pp_ect.load_CDFfiles_ECT(sdate_ect, edate_ect, data_dir,
                     ect_var, ect_rename, probe_opt, 'mageis', '3', 'fedu')

    [rept_info, rept_fedu, rept_info_meta, rept_fedu_meta] = rept_probe
    [mageis_info, mageis_fedu, mageis_info_meta, mageis_fedu_meta] = mageis_probe

    ### Promediamos y sincronizamos
    rept_info, rept_fedus = mf.time_average_ect_data(rept_info, rept_fedu,
                            sdate_ect, edate_ect, time_avg)
    mageis_info, mageis_fedus = mf.time_average_ect_data(mageis_info,
                                mageis_fedu, sdate_ect, edate_ect, time_avg)

    ect_info = mf.sync_ect_data(rept_info, mageis_info)

    ### Unidades para energy and alpha bins for REPT and MagEIS
    rept_energy_bins =  rept_fedu_meta['energy_values']*u.MeV
    rept_alpha_bins = rept_fedu_meta['alpha_values']*u.degree
    rept_bins = [rept_energy_bins, rept_alpha_bins]

    mageis_energy_bins =  mageis_fedu_meta['energy_values']*u.keV
    mageis_alpha_bins = mageis_fedu_meta['alpha_values']*u.degree
    mageis_bins = [mageis_energy_bins, mageis_alpha_bins]

    return ect_info, [rept_fedus, rept_bins], [mageis_fedus, mageis_bins]


def step0_get_OMNI_data(stime, etime, data_dir, resolution, file_type,
                        flag_download, flag_model):
    print('------------------------------------------------------------------------')
    print('STEP 0: Get and process OMNI data')
    print(f'    Only return relevant keys for magnetic field model {flag_model}')
    print('------------------------------------------------------------------------')

    ### Fechas descarga ect
    sdate_ect = stime[:10]
    edate_ect = etime[:10]

    ### Fechas descarga omni
    sdate_omni_1h, sdate_omni_min, edate_omni = mf.get_omni_dates(sdate_ect, edate_ect)

    if flag_download:
        remote_dir = "https://cdaweb.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/"

        # OMNI
        dd_omni.download_CDFfiles_OMNI(sdate_omni_1h, edate_omni, remote_dir,
                    data_dir, '1h', '')

        if resolution !='1h':
            dd_omni.download_CDFfiles_OMNI(sdate_omni_min, edate_omni, remote_dir,
                    data_dir, resolution, file_type)
    else:
        print('* No download')

    ### Variables omni según resolución y modelo de campo magnético
    var_1h, var_min, rename, range = mf.model_variables_omni(flag_model, resolution)

    ### OMNI
    omni_1h = pp_omni.load_CDFfiles_OMNI(sdate_omni_1h, edate_omni, data_dir,
                      var_1h, rename, '1h', '')

    omni_min = [False, False]
    if resolution !='1h':
         omni_min = pp_omni.load_CDFfiles_OMNI(sdate_omni_min, edate_omni,
                            data_dir, var_min, rename, resolution, file_type)

    ### Juntamos info omni resolución de 1 hour & min
    omni_info, omni_meta = mf.sync_omni_data(omni_1h, omni_min, resolution)

    return omni_info, omni_meta, range


def step0_inputs(ect_info, omni_info, stime, etime, res_omni, range_omni,
    rept_fedus, mageis_fedus, units):

    units_inv, units_flux, _ = units
    Re, _, _, _ = units_inv
    unit_flux_mev, unit_flux_kev = units_flux

    # Datos ect
    ect_filt = mf.filter_dates(ect_info, stime, etime)
    idx1 = ect_filt.index[0]
    idx2 = ect_filt.index[-1]
    dates, xs_geo, ys_geo, zs_geo = mf.get_GEO_coordinates(ect_filt, Re)
    list_x_inputs, N = mf.dicts_x_input(dates, xs_geo, ys_geo, zs_geo)

    # Datos omni
    delta_omni = pd.Timedelta(res_omni)
    omni_filt = mf.filter_dates(omni_info, stime, etime, delta_omni)
    omni_filt = mf.filter_mf_inputs(omni_filt, range_omni)

    list_mag_inputs = mf.dicts_magnetic_input(dates, omni_filt)

    # Unidades para fedus de REPT y MagEIS
    rept_fedus_inputs = rept_fedus[idx1:idx2+1]*unit_flux_mev
    mageis_fedus_inputs = mageis_fedus[idx1:idx2+1]*unit_flux_kev

    return [list_x_inputs, list_mag_inputs], [rept_fedus_inputs, mageis_fedus_inputs], N


def get_inputs(times, sat, data_dirs, flag_dwnl, opts_model, time_avg, res_omni,
               file_type_omni, targets, units):

    start_time, end_time = times
    data_dir_ect, data_dir_omni = data_dirs
    options_mf, flag_model = opts_model
    target_Ks, target_mus = targets


    ### Descarga y lectura de datos ECT-REPT, ECT-MagEIS, y OMNI
    # ECT
    ect_info, rept, mageis = step0_get_ECT_data(start_time, end_time, sat,
                             data_dir_ect, flag_dwnl, flag_model, time_avg)

    rept_fedus, rept_bins = rept
    mageis_fedus, mageis_bins = mageis

    # OMNI
    omni_info, omni_meta, omni_range = step0_get_OMNI_data(start_time, end_time,
                           data_dir_omni, res_omni, file_type_omni, flag_dwnl,
                           flag_model)

    ### Creamos los inputs
    inputs = step0_inputs(ect_info, omni_info, start_time, end_time, res_omni,
             omni_range, rept_fedus, mageis_fedus, units)

    return inputs, [rept_bins, mageis_bins]



def step1(mf_model, inputs, alphas, info_K, unit_K):
    # Step 1: Calculamos K para distintos alphas
    print('------------------------------------------------------------------------')
    print('STEP 1: K calculation from IRBEM for different PA values')
    print('------------------------------------------------------------------------')

    func_K, str_func, flag_K = info_K

    print(f'* Using option {flag_K} for calculation:')
    print(f'      - {str_func}')

    # Obtenemos los valores de K para el arreglo de PA en unidades de Re*(nT**(1/2))
    with open(os.devnull, 'w') as fnull:
        with contextlib.redirect_stdout(fnull), contextlib.redirect_stderr(fnull):
            Ks, _ = func_K(alphas, mf_model, inputs, unit_K)
#    Ks, _ = func_K(alphas, mf_model, inputs, unit_K)
    ### Limpiamos datos nan
    mask = np.isnan(Ks.value)
    Ks = Ks[~mask]
    alphas = alphas[~mask]
    print()
    return alphas, Ks


def step2(alphas, Ks, target_Ks):
    # Step 2: Interpolamos la función alpha(K) para obtener target alphaK
    print('------------------------------------------------------------------------')
    print('STEP 2: Interpolate function alpha(K) to get target alphaK for fixed K')
    print('------------------------------------------------------------------------')

    # Obtenemos el interpolador spline para alpha(K)
    spline = inv.interpolator_alpha(alphas, Ks)

    # Entre estos valores puedo elegir K
    K_min = Ks[-1]
    K_max = Ks[0]

    print('* K range for this time step:')
    print(f'      - [{K_min}, {K_max}]')

    # Obtenemos el valor de alpha para cada uno de los target K
    target_alphasK = np.empty(len(target_Ks))
    for i, K in enumerate(target_Ks):
        alphaK = inv.interpolate_alpha_K(K, K_max, spline)
        target_alphasK[i] = alphaK.value

    out_alphas = target_alphasK*alphas.unit
    print('* Target Ks: ', target_Ks)
    print('* Target alphasK: ', out_alphas)
    print()
    return out_alphas, spline, K_min, K_max


def step3_model(mf_model, inputs, alphasK, mus):
    # Step 3: Calcular target energy(mu,K) of chosen mu and K
    print('------------------------------------------------------------------------')
    print('STEP 3: Calculate target E_mu for fixed mu and K')
    print('------------------------------------------------------------------------')
    print('* Using model for local magnenitc field for calculation:')
    print()
    list_Esmu = []
    ### Obtenemos el campo magnético local en la posición del spacecraft
    ### usando IRBEM, recordar que está en nanoteslas
    with open(os.devnull, 'w') as fnull:
        with contextlib.redirect_stdout(fnull):
            dict_mag_field = mf_model.get_field_multi(*inputs)
            b_mag = dict_mag_field['Bl'][0]*u.nT
#    print(b_mag)
    for alpha in alphasK:
        Esmu = inv.calculate_E(mus, b_mag, alpha)
        list_Esmu.append(Esmu)
        print('* Target mus: ', mus)
        print('* Target alphaK: ', alpha)
        print('* Target Esmu: ', Esmu)
        print()
    return list_Esmu


def step4(mf_model, inputs, alphasK, info_Lstar, save_lstar):
    # Step 4: Calculamos Lstar para partículas que tienen PA alphaK
    print('------------------------------------------------------------------------')
    print('STEP 4: Lstar calculation from IRBEM for alphaK')
    print('------------------------------------------------------------------------')

    func_Lstar, str_func, flag_Lstar = info_Lstar
    df_lstar, index, Ks, mus = save_lstar

    print(f'* Using option {flag_Lstar} for calculation:')
    print(f'      - {str_func}')

    N_mus = len(mus)
    with open(os.devnull, 'w') as fnull:
        with contextlib.redirect_stdout(fnull):
            for i, alpha in enumerate(alphasK):
                if np.isnan(alpha):
                    lstar=np.nan
                else:
                    lstar = func_Lstar(alpha, mf_model, inputs)
                    if lstar <0:
                        lstar = np.nan
                df_lstar.loc[index, str(Ks[i])] = [lstar]*N_mus

    print('* Target alphasK: ', alphasK)
    print('* Lstar at target alphasK:\n ', df_lstar.loc[index])
    print()
    return df_lstar


def step5_onefunc(fit_info, fedu, channels_to_use, list_bins, N_energy, alphasK, inst, save_errors):
    # Step 5: Interpolar datos j(alpha) para cada canal de energía para obtener
    # flujos en target alphaK
    print('-------------------------------------------------------------------')
    print('STEP 5: Interpolate j(alpha). Get flux at target alphaK for each')
    print('                       energy channel')
    print('-------------------------------------------------------------------')

    func_opt, str_func, func, p0, lims = fit_info
    energy_bins, alpha_bins = list_bins
    df_rmse, df_r2, index = save_errors
    cols = df_rmse.columns[df_rmse.columns.str.contains(inst)]

    print(f'* Using function {func_opt} for interpolation:')
    print(f'      {str_func}')
    print(f'* Instrumet: {inst}')
    # Obtenemos el mejor ajuste

    if np.all(np.isnan(alphasK)):
#        print('Todos los alphas son nan, así que no vale la pena hacer los ajustes')
        js_alphasK = np.array([[]]*len(alphasK)).T

        df_rmse.loc[index, cols] = np.full(len(channels_to_use), np.nan)
        df_r2.loc[index, cols] = np.full(len(channels_to_use), np.nan)

        return js_alphasK*fedu.unit, []*energy_bins.unit, 0

    else:
#        print('Al menos un alpha no es nan, así que tenemos que hacer los ajustes')
        with open(os.devnull, 'w') as fnull:
            with contextlib.redirect_stdout(fnull):
                fit_results = fp.fitting_alpha_flux(func, p0, lims, fedu, channels_to_use, alpha_bins, N_energy)
        fit, parameters, max_values= fit_results[:3]
        err_fit = fit_results[4]

        df_rmse.loc[index, cols] = err_fit[:,0]
        df_r2.loc[index, cols] = err_fit[:,1]

        # Usamos el mejor ajuste para calcular flujo at target alphaK, j(alphaK)
        array_indexes = np.array(range(N_energy))[fit]
        js_alphasK = np.full((len(array_indexes), len(alphasK)), np.nan)
        for i, idx in enumerate(array_indexes):
            max_flux = max_values[idx]
            par_opt = parameters[idx]
            js_idx = func(np.radians(alphasK.value), *par_opt)
            js_alphasK[i, :] = js_idx*max_flux.value

        energy_bins_to_use = energy_bins[fit]
        print()
        return js_alphasK*max_values.unit, energy_bins_to_use, fit_results



def step5_bestfunc(list_fit_info, fedu, channels_to_use, list_bins, N_energy, alphasK, inst, save_errors):
    # Step 5: Interpolar datos j(alpha) para cada canal de energía para obtener
    # flujos en target alphaK
    print('-------------------------------------------------------------------')
    print('STEP 5: Interpolate j(alpha). Get flux at target alphaK for each')
    print('                       energy channel')
    print('-------------------------------------------------------------------')

    energy_bins, alpha_bins = list_bins
    df_rmse, df_r2, index = save_errors
    cols = df_rmse.columns[df_rmse.columns.str.contains(inst)]

    print(f'* Using functions best function for interpolation:')
    print(f'      ')
    print(f'* Instrumet: {inst}')

    if np.all(np.isnan(alphasK)):
#        print('Todos los alphas son nan, así que no vale la pena hacer los ajustes')
        js_alphasK = np.array([[]]*len(alphasK)).T

        df_rmse.loc[index, cols] = np.full(len(channels_to_use), np.nan)
        df_r2.loc[index, cols] = np.full(len(channels_to_use), np.nan)

        return js_alphasK*fedu.unit, []*energy_bins.unit, 0

    else:
#        print('Al menos un alpha no es nan, así que tenemos que hacer los ajustes')

        # Obtenemos el mejor ajuste
        with open(os.devnull, 'w') as fnull:
            with contextlib.redirect_stdout(fnull):
                fit_results = fp.fitting_alpha_flux_V2(list_fit_info, fedu, channels_to_use, alpha_bins, N_energy)

        fit, parameters, max_values, _, err_fit, fit_opts, funcs = fit_results

        df_rmse.loc[index, cols] = err_fit[:,0]
        df_r2.loc[index, cols] = err_fit[:,1]

        # Usamos el mejor ajuste para calcular flujo at target alphaK, j(alphaK)
        array_indexes = np.array(range(N_energy))[fit]
        js_alphasK = np.full((len(array_indexes), len(alphasK)), np.nan)
        for i, idx in enumerate(array_indexes):
            func = funcs[idx]
            max_flux = max_values[idx]
            par_opt = parameters[idx]
            js_idx = func(np.radians(alphasK.value), *par_opt)
            js_alphasK[i, :] = js_idx*max_flux.value

        energy_bins_to_use = energy_bins[fit]
        print()
        return js_alphasK*max_values.unit, energy_bins_to_use, fit_results



def step6_flux(energy, array_flux, fit_info, Esmu, alphasK, units_psd, save_psd):
    # Step 6: Obtain f(Emu) using calculated j(alphaK) for each energy channel,
    # to interpolate j(E)
    print('-------------------------------------------------------------------')
    print('STEP 6: Interpolate j(E) using calculated fluxes at target alphasK.')
    print('            Get psd at target Emu for each alphaK')
    print('-------------------------------------------------------------------')

    fit_opt, str_func, fitter, p0, lims = fit_info

    df_psd, index, Ks = save_psd
    print(f'* Using function {fit_opt} for interpolating flux:')
    print(f'      {str_func}')
    try:
        E_max = energy[-1]
    except:
        E_max = 0

#    print('E_max', E_max)
#    print('target Esmu', Esmu)
    list_params = []
    list_fit_obj = []
    for i in range(len(alphasK)):
        flux = array_flux[:,i]
        Es = Esmu[i]
        mask = np.isnan(flux.value)
        flux = flux[~mask]
#        print(i, flux, Es)
        Es2 = Es.copy()
        if E_max ==0:
            # No hay fit
            j_Es = np.array([np.nan]*len(Es))*flux.unit
            fit_obj = False
            parms = False
        elif np.all(np.isnan(Es)):
            # No hay fit
            j_Es = np.array([np.nan]*len(Es))*flux.unit
            fit_obj = False
            parms = False

        else:
            bool2 = Es2 > E_max
#            print(bool2)
            Es2[bool2] = np.nan*E_max.unit
#            print(Es2)

            with open(os.devnull, 'w') as fnull:
                with contextlib.redirect_stdout(fnull):
                    fit_results = fp.fitted_y_at_Emu(energy, flux, fitter, p0, lims,'flux', Es2)
            fit_obj, parms, j_Es = fit_results

        f_Es = fp.flux_to_psd(Es, j_Es, units_psd)
#        print(j_Es, f_Es)
        df_psd.loc[index, str(Ks[i])] = f_Es
        list_params.append(parms)
        list_fit_obj.append(fit_obj)

        print('* Target K: ', Ks[i])
        print('* Target alphaK: ', alphasK[i])
        print('* Target Esmu: ', Es2)
        print('* PSD at target alphaK and Esmu: ', f_Es)
        print()
    return df_psd, list_fit_obj, list_params


def step6_psd(energy, array_flux, fit_info, Esmu, alphasK, units_psd, save_psd):
    # Step 6: Obtain f(Emu) using calculated j(alphaK) for each energy channel,
    # to interpolate f(E)
    print('-------------------------------------------------------------------')
    print('STEP 6: Interpolate f(E) using calculated fluxes at target alphasK.')
    print('            Get psd at target Emu for each alphaK')
    print('-------------------------------------------------------------------')

    fit_opt, str_func, fitter, p0, lims = fit_info
    df_psd, index, Ks = save_psd
    print(f'* Using function {fit_opt} to interpolate psd:')
    print(f'      {str_func}')

    ### Calculamos psd to fit at target alphaK
    try:
        E_max = energy[-1]
    except:
        E_max = 0
#    print('E_max', E_max)
#    print('target Esmu', Esmu)

    array_psd = fp.flux_to_psd(energy[:,np.newaxis], array_flux, units_psd)
#    print(array_psd)
    list_params = []
    list_fit_obj = []
    for i in range(len(alphasK)):
        psd = array_psd[:,i]
        Es = Esmu[i]
        mask = np.isnan(psd.value)
        psd = psd[~mask]
#        print(i, psd, Es)
        Es2 = Es.copy()

        if E_max ==0:
            # No hay fit
            f_Es = np.array([np.nan]*len(Es))*psd.unit
            fit_obj = False
            parms = False

        elif np.all(np.isnan(Es)):
            # No hay fit
            f_Es = np.array([np.nan]*len(Es))*psd.unit
            fit_obj = False
            parms = False

        else:
            bool2 = Es2 > E_max
#            print(bool2)
            Es2[bool2] = np.nan*E_max.unit
#            print(Es2)


            with open(os.devnull, 'w') as fnull:
                with contextlib.redirect_stdout(fnull):
                    fit_results = fp.fitted_y_at_Emu(energy, psd, fitter, p0, lims,'psd', Es2)
            fit_obj, parms, f_Es = fit_results

        df_psd.loc[index, str(Ks[i])] = f_Es
        list_params.append(parms)
        list_fit_obj.append(fit_obj)

        print('* Target K: ', Ks[i])
        print('* Target alphaK: ', alphasK[i])
        print('* Target Esmu: ', Es2)
        print('* PSD at target alphaK and Esmu: ', f_Es)
        print()

    return df_psd, list_fit_obj, list_params


def psd_calculation(channels_to_use, options_psd, options_model, targets,
                    inputs, bins, units, N_steps = 10e6, N_alphas = 20):

    rept_cut_offs, rept_N, mageis_cut_offs, mageis_N = channels_to_use
    PA_fit_opt, energy_fit_opt, K_opt, Lstar_opt = options_psd

    options_mf, flag_model = options_model
    target_Ks, target_mus = targets
    inputs_mf, inputs_fedus, N_data = inputs
    rept_bins, mageis_bins = bins
    units_inv, units_flux, units_psd = units

    # Array de PAs en grados
    alphas_right = np.geomspace(30, 90, N_alphas +1)[1:]*u.degree
    alphas_left = np.geomspace(30, 3, N_alphas)[::-1]*u.degree
    alphas_total = np.concatenate([alphas_left, alphas_right])

    _, _, unit_K, _ = units_inv
    unit_flux_mev, _ = units_flux

    ### Creamos la lista de los canales de energía a usar
    rept_channels = fp.channels_to_use(rept_cut_offs, rept_N)
    mageis_channels = fp.channels_to_use(mageis_cut_offs, mageis_N)

    ### Obtenemos info para los fits de pitch angle.
    PA_fit_info = fp.info_fit_PA_flux(PA_fit_opt)
    step5 = step5_onefunc
    if PA_fit_info[2] == None:
        PA_fit_info = [fp.info_fit_PA_flux(x) for x in PA_fit_opt]
        step5 = step5_bestfunc
    print(PA_fit_info)
    print(step5)

    ### Obtenemos info para los fits de energía.
    energy_fit_info = fp.info_fit_energy_data(energy_fit_opt)

    ### Obtenemos función para el cálculo de K y Lstar
    info_K = inv.info_calculate_K(K_opt)
    info_Lstar = inv.info_calculate_Lstar(Lstar_opt)

    ### Fijamos cuál función para el step6 vamos a usar
    # Depende de si queremos interpolar f(E) o j(E)
    if energy_fit_opt[1] == 'flux':
        step6 = step6_flux
    elif energy_fit_opt[1] == 'psd':
        step6 = step6_psd

    ### Creamos los dataframes para guardar datos
    labels_Ks = [str(x) for x in target_Ks]
    labels_mus = [str(x) for x in target_mus]

    cols = pd.MultiIndex.from_product([labels_Ks, labels_mus], names=['K', 'mu'])
    df_psd = pd.DataFrame(columns = cols)
    df_lstar = pd.DataFrame(columns=cols)

    ### Creamos los dataframes para guardar errores
    rept_cols = [f'REPT_{channel}' for channel in range(1, rept_N+1)]
    mageis_cols = [f'MagEIS_{channel}' for channel in range(1, mageis_N+1)]
    all_cols = rept_cols + mageis_cols
    df_rmse = pd.DataFrame(columns=all_cols)
    df_r2 = pd.DataFrame(columns=all_cols)

    ### Obtenemos los inputs para el cálculo de psd
    x_inputs, mag_inputs = inputs_mf
    rept_fedus_inputs, mageis_fedus_inputs = inputs_fedus

    ### Creamos el objeto de campo magnético
    model_obj = IRBEM.MagFields(options=options_mf, kext=flag_model,
                                verbose = False, sysaxes=1)
    N_steps = min(N_steps, N_data)
    for i in range(N_steps):
        # Obtenemos los inputs del modelo de campo magnético, para el timestep i
        inputs = [x_inputs[i], mag_inputs[i]]
        timestamp = x_inputs[i]['dateTime']
        print(f'TIMESTAMP {i}: ', timestamp)

        # Obtenemos los flujo de REPT y MagEIS para el timestep i
        rept_fedu = rept_fedus_inputs[i]
        mageis_fedu = mageis_fedus_inputs[i]

        save_psd = [df_psd, timestamp, target_Ks]
        save_lstar = [df_lstar, timestamp, target_Ks, target_mus]
        save_errors = [df_rmse, df_r2, timestamp]

        ###################  Section 1: Calculo de invariantes  ################

        ''' Step 1: Calcular K para distintos PA alphas '''
        alphas, Ks = step1(model_obj, inputs, alphas_total, info_K, unit_K)
#        print(alphas)
#        print(Ks)
        ''' Step 2: Interpolar la función alpha(K) y calcular target alpha_K '''
        target_alphasK, spline, K_min, K_max = step2(alphas, Ks, target_Ks)

        ''' Step 3: Calcular energía(mu,K) of chosen mu and K '''
        target_Esmu = step3_model(model_obj, inputs, target_alphasK, target_mus)

        ''' Step 4: Calculate L* '''
        df_lstar = step4(model_obj, inputs, target_alphasK, info_Lstar,
                   save_lstar)


        #######################  Section 2: FEDU processing  ###################

        ''' Step 5: Obtain j(alphaK) for each energy channel '''
        ### Ajustar REPT y MagEIS PA flux y calculamos flux at target alphaK
        rept = step5(PA_fit_info, rept_fedu, rept_channels, rept_bins, rept_N,
               target_alphasK, 'REPT', save_errors)
        mageis  = step5(PA_fit_info, mageis_fedu, mageis_channels, mageis_bins,
                  mageis_N, target_alphasK, 'MagEIS', save_errors)

        rept_flux_alphasK, rept_energy, rept_PA_fit_res = rept
        mageis_flux_alphasK, mageis_energy, mageis_PA_fit_res = mageis


        flux_alphasK = [rept_flux_alphasK, mageis_flux_alphasK]
        energy_bins = [rept_energy, mageis_energy]

        ''' Step 6: Obtain psd(Emu) using j(alphaK) for each energy '''
        ### Mezclamos mageis y rept fluxes at target alphaK
        energy_to_fit, energy_flux_to_fit = fp.join_energy_flux(flux_alphasK,
                                            energy_bins, unit_flux_mev)

        ### Ajustar y=flux/psd at target alphaK y calculamos psd at target E_mu
        df_psd, _, _ = step6(energy_to_fit, energy_flux_to_fit, energy_fit_info,
                       target_Esmu, target_alphasK, units_psd, save_psd)

    return df_psd, df_lstar, df_rmse, df_r2
