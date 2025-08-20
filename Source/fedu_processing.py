#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import curve_fit
from astropy import units as u
from astropy import constants as const
import scipy.interpolate as scipy_int


'''
###############################################################################
### En este script están todas las funciones para hacer el procesamiento de
### flujos, necesario para el cálculo PSD
################################################################################
'''

#################################### PA fit ####################################

def channels_to_use(cut_offs, N_energy):
    channels = [False]*(cut_offs[0]-1) + [True]*(cut_offs[1] -
               cut_offs[0]+1) + [False]*(N_energy-cut_offs[1])
    return channels


def PA_flux_func1(alpha_rad, c0, c1, n):
    j = c0*np.sin(alpha_rad) + c1*(np.sin(alpha_rad)**n)
    return j


def PA_flux_func2(alpha_rad, c, n):
    j = c*(np.sin(alpha_rad)**n)
    return j


def info_fit_PA_flux(flag_fit):
    if flag_fit =='1':
        func = PA_flux_func1
        str_func = '$j(\\alpha) = c_0 sin(\\alpha) + c_1 sin^n(\\alpha)^n$'
        p0 = [1,0.1,1]
        lims = ([0,-2,0], [2, 2, 100])

    elif flag_fit =='2':
        func = PA_flux_func2
        str_func = '$j(\\alpha) = C sin^n(\\alpha)$'
        p0 = [2,1]
        lims = ([0,0], [10, 100])
    else:
         func = None
         str_func = ''
         p0 = None
         lims = None
    print('Fuction %s to fit PA flux\n%s' %(flag_fit, str_func))

    return flag_fit, str_func, func, p0, lims


def fitting_alpha_flux(func, p0, lims, fedu, channels_to_use, alpha_bins, N_energy):
    fit = []
    parameters = []
    max_values = []
    err_par = []
    err_fit = []
    fit_opts = []
    for i in range(1, N_energy+1):
        print('------')
        print('Energy channel %s' % i)
        if not channels_to_use[i-1]:
            print('No vamos a usar este canal de energía')
            fit.append(False)
            parameters.append(np.nan)
            err_par.append(np.nan)
            err_fit.append(np.nan)
            max_values.append(np.nan)
        else:
            print('Si vamos a usar este canal de energía')
            flux = fedu[:,i-1]
            # Eliminamos los nan del flux
            mask_flux = np.isnan(flux)
            flux_red = flux[~mask_flux]
            alpha_bins_red = alpha_bins[~mask_flux]
            try:
                max_flux = np.max(flux_red)
                max_values.append(max_flux.value)
            except:
                max_flux = np.nan
                max_values.append(max_flux)

            if max == 0.:
                print('No fit: fluxes are all 0')
                fit.append(False)
                parameters.append(0)
                err_par.append(np.nan)
                err_fit.append(np.nan)

                continue
            elif np.isnan(max_flux):
                print('No fit: fluxes are all nan')
                fit.append(False)
                parameters.append(0)
                err_par.append(np.nan)
                err_fit.append(np.nan)

                continue
            else:
                flux_norm = flux_red/max_flux
            # Hacemos el ajuste solo si hay 4 o más valores distintos de 0
            if len(np.nonzero(flux_red)[0]) > 3:
                par_opt, pcov = curve_fit(func, np.radians(alpha_bins_red), flux_norm, p0=p0, bounds = lims)
                print('Fit: At least 4 usable flux values')
                fit.append(True)
                parameters.append(par_opt)
                err_par.append(np.sqrt(np.diag(pcov)))
                err_fit.append(errors(np.radians(alpha_bins_red), flux_norm, par_opt, func))
            else:
                print('No fit: Less than 4 usable flux values')
                fit.append(False)
                parameters.append(0)
                err_par.append(np.nan)
                err_fit.append(np.nan)

    return np.array(fit), parameters, max_values*fedu.unit, err_par, err_fit


def fitting_alpha_flux_V2(list_info_fit, fedu, channels_to_use, alpha_bins, N_energy):

    fit = []
    parameters = []
    max_values = []
    err_par = []
    err_fit = []
    fit_opts = []
    fit_func = []

    for i in range(1, N_energy+1):
        print('------')
        print('Energy channel %s' % i)
        if not channels_to_use[i-1]:
            print('No vamos a usar este canal de energía')
            fit.append(False)
            parameters.append(np.nan)
            max_values.append(np.nan)
            err_par.append(np.nan)
            err_fit.append(np.nan)
            fit_opts.append(False)
            fit_func.append(False)
        else:
            print('Si vamos a usar este canal de energía')
            flux = fedu[:,i-1]
            # Eliminamos los nan del flux
            mask_flux = np.isnan(flux)
            flux_red = flux[~mask_flux]
            alpha_bins_red = alpha_bins[~mask_flux]
            try:
                max_flux = np.max(flux_red)
                max_values.append(max_flux.value)
            except:
                max_flux = np.nan
                max_values.append(max_flux)

            if max == 0.:
                print('No fit: fluxes are all 0')
                fit.append(False)
                parameters.append(0)
                err_par.append(np.nan)
                err_fit.append(np.nan)
                fit_opts.append(False)
                fit_func.append(False)
                continue
            elif np.isnan(max_flux):
                print('No fit: fluxes are all nan')
                fit.append(False)
                parameters.append(0)
                err_par.append(np.nan)
                err_fit.append(np.nan)
                fit_opts.append(False)
                fit_func.append(False)
                continue
            else:
                flux_norm = flux_red/max_flux
            # Hacemos el ajuste solo si hay 4 o más valores distintos de 0
            if len(np.nonzero(flux_red)[0]) > 3:
                print('Fit: At least 4 usable flux values')
                fit.append(True)
                res_error = []
                for info_fit in list_info_fit:
                    opt, str_title, func, p0, lims = info_fit
                    par_opt, pcov = curve_fit(func, np.radians(alpha_bins_red), flux_norm, p0=p0, bounds = lims)
                    rmse, r2 = errors(np.radians(alpha_bins_red), flux_norm, par_opt, func)
                    res_error.append({'r2':r2.value, 'rmse':rmse.value, 'opt':opt, 'par_opt':par_opt, 'func':func})

                best = max(res_error, key=lambda x: x['r2'])
                parameters.append(best['par_opt'])
                err_par.append('miau')
#                err_par.append(np.sqrt(np.diag(pcov)))
                err_fit.append([best['r2'], best['rmse']])
                fit_opts.append(best['opt'])
                fit_func.append(best['func'])
            else:
                print('No fit: Less than 4 usable flux values')
                fit.append(False)
                parameters.append(0)
                err_par.append(np.nan)
                err_fit.append(np.nan)
                fit_opts.append(False)
                fit_func.append(False)

    return np.array(fit), parameters, max_values*fedu.unit, err_par, err_fit, fit_opts, fit_func


def errors(x,y, par_opt, func):
    N = len(y)
    p = len(par_opt)
    res = y - func(x, *par_opt)
    rss = np.sum(res**2)
    tss = np.sum((y-np.mean(y))**2)
    rmse = np.sqrt(rss/N)
    r2 = 1.- rss/tss
    return rmse, r2


def fitted_flux_at_alphaK(fit_results, list_func, alphaK):
    fit, parameters, max_values = fit_results[:3]
    N = len(fit)
    array_indexes = np.array(range(N))[fit]
    js_alphaK = np.full(N, np.nan)
    for idx in array_indexes:
        func = list_func[idx]
        max_flux = max_values[idx]
        par_opt = parameters[idx]
        j_idx = func(np.radians(alphaK.value), *par_opt)
        js_alphaK[idx] = j_idx*max_flux.value

    js_alphaK = js_alphaK[~np.isnan(js_alphaK)]
    return js_alphaK*max_values.unit


################################## Energy fit ##################################

def join_energy_flux(list_flux_to_fit, list_bins, unit_flux_mev):

    rept_flux_to_fit, mageis_flux_to_fit = list_flux_to_fit
    rept_energy_to_fit, mageis_energy_to_fit = list_bins

    rept_flux = rept_flux_to_fit.to(unit_flux_mev)
    mageis_flux = mageis_flux_to_fit.to(unit_flux_mev)
    flux = np.concatenate((rept_flux, mageis_flux), axis=0)

    rept_energy = rept_energy_to_fit.to(u.MeV)
    mageis_energy = mageis_energy_to_fit.to(u.MeV)
    energy = np.concatenate((rept_energy, mageis_energy), axis=0)

    sorted_index = np.argsort(energy.value)
    sorted_energy = energy.value[sorted_index]*energy.unit
    sorted_flux = flux.value[sorted_index]*flux.unit

    return sorted_energy, sorted_flux


def flux_to_psd(energy, j, units_psd):
    unit_c, unit_psd = units_psd
#    flux = j.to(unit_flux_mev)
    flux = j
#    print('flux', flux)
    E= energy.to(u.MeV)
#    print('E', E)
    c2m0 = (const.m_e*const.c**2).to('MeV')
#    print('c2m0', c2m0)
    p2 = (E**2 + 2*c2m0*E)/(unit_c**2)
#    p2 = p2[:, np.newaxis]
#    print('p2', p2)
    flux_c  = (flux*u.sr)/(const.c.to('cm/s'))
#    print('flux_c', flux_c)
    f = (flux_c/p2)*unit_c

    psd = f.to(unit_psd)
#    print('f', f)
    return psd


def energy_flux_pow_law(log_energy, c0, c1):
    log_j = c0*log_energy + c1
    return log_j


def energy_psd_exp(energy, c0, c1):
    f = c0*np.exp(-c1*energy)
    return f


def fitter_energy_flux_pl(data_type, energy, flux, p0, lims, Emu):
    print('\nFitting using flux, pl')
    func  = energy_flux_pow_law
    par_opt, _ = curve_fit(func, np.log10(energy.value), np.log10(flux.value), p0=p0, bounds = lims)
    log_j_Emu = func(np.log10(Emu.value), *par_opt)
    j_Emu = (10**log_j_Emu)*flux.unit
    out = func
    return out, par_opt, j_Emu


def fitter_energy_psd_exp(data_type, energy, psd, p0, lims, Emu):
    print('\nFitting using psd, exp')
    func = energy_psd_exp
    par_opt, _ = curve_fit(func, energy.value, psd.value, p0=p0, bounds = lims)
    f_Emu = func(Emu.value, *par_opt)
    out = func
    return out, par_opt, f_Emu*psd.unit


def fitter_energy_y_spl(data_type, energy, y, p0, lims, Emu):
    print('\nFitting using %s, spl' %data_type)

    w = np.ones_like(energy.value)
#    w[:6] += 5.5
#    w[6:] -= 0.5
    print(w)
    spl_obj = spy_int.make_smoothing_spline(np.log10(energy.value), np.log10(y.value), w=w)
#    spl_obj = spy_int.UnivariateSpline(np.log10(energy.value), np.log10(data.value), s=0, w=w)
    log_y_Emu = spl_obj(np.log10(Emu.value))
    y_Emu = 10**log_y_Emu*y.unit
    out = spl_obj
    par_opt = []
    return out, par_opt, y_Emu


def fitter_energy_y_lin(data_type, energy, y, p0, lims, Emu):
    print('\nFitting using %s, lin' %data_type)
    spl_lin_obj = scipy_int.make_interp_spline(np.log10(energy.value), np.log10(y.value), k=1)
    log_y_Emu = spl_lin_obj(np.log10(Emu.value))
    y_Emu = 10**log_y_Emu*y.unit
    out = spl_lin_obj
    par_opt = []
    return out, par_opt, y_Emu


def info_fit_energy_data(fit_opt):
    if fit_opt[0] == 'pl':
        fitter = fitter_energy_flux_pl
        str_func = 'log(j(E)) = c0*log(E)+ c1'
        p0 = [-10,1]
        lims = ([-np.inf,0], [0, np.inf])
    elif fit_opt[0] == 'spl':
        fitter = fitter_energy_y_spl
        str_func = 'Smoothing spline'
        p0 = False
        lims = False
    elif fit_opt[0] == 'exp':
        fitter = fitter_energy_psd_exp
        str_func = 'f(E) = c0*exp**(-c1*E)'
        p0 = [10,1]
        lims = ([0,0], [np.inf, np.inf])
    elif fit_opt[0] == 'lin':
        fitter = fitter_energy_y_lin
        str_func = 'Linear interpolator'
        p0 = False
        lims = False
    print('Fuction %s to fit energy %s data\n%s' %(fit_opt[0], fit_opt[1], str_func))
    return fit_opt, str_func, fitter, p0, lims


def fitted_y_at_Emu(energy, y, fitter, p0, lims, data_type, Esmu):
    N = len(Esmu)
    if len(y)< 6:
        print('------\nLess than 6 points to fit')
        y_Esmu = np.array([np.nan]*N)*y.unit
        fit_obj = False
        parms = False
    else:
        print('------\n6 or more points to fit')
        fit_obj, parms, y_Esmu = fitter(data_type, energy, y, p0, lims, Esmu)

    return fit_obj, parms, y_Esmu
