#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from astropy import units as u
import os
import contextlib
import matplotlib.dates as mdates

import Source.plots_fedu_processing as plots_fp
import Source.plots_invariants_calculation as plots_inv

#matplotlib.rcParams['text.usetex'] = True
#plt.rcParams['font.family'] = 'serif'

'''
###############################################################################
### En este script están todas las funciones para visualizar que los cálculos
### que involucran invariantes adiabáticas sean correctos
################################################################################
'''


def check_step2(info_K, model_obj, inputs, unit_K, Ks_spline, alphas_spline, spline,
                Ks, alphasK):
    func_K, _, _ = info_K

    alphas_plot = np.linspace(0, 90, 100)*u.degree
    with open(os.devnull, 'w') as fnull:
        with contextlib.redirect_stdout(fnull):
            Ks_plot, _ = func_K(alphas_plot, model_obj, inputs, unit_K)
    list_Ks = [Ks_plot, Ks_spline]
    list_alphas = [alphas_plot, alphas_spline]
    list_labels = ['Calculated with MF model', 'Used for spline']
    list_colors = ['silver', 'b']
    ax = plots_inv.plot_alpha_K(list_Ks, list_alphas, list_labels, list_colors)
    ax.plot(np.log10(Ks.value), alphasK.value, 'o', color = 'r', label = 'Target alphasK')
    ax.plot(np.log10(Ks_plot.value), spline(np.log10(Ks_plot.value)), color = 'r')
    ax.legend()
    plt.show()


def check_step5(fit_info, fedu, energy_bins, fit_results, js_alphasK, alphasK, instrument):

    plots_fp.check_fit_PA_flux(fit_info, fedu, energy_bins, fit_results,
                                js_alphasK, alphasK, instrument, '', 1, 0)


def check_step6(js, energy, alphasK, Ks, fit_objects, parameters, df_psd, fit_info,
                energy_range, Esmu, unit_c, unit_flux, func):

    for i, alphaK in enumerate(alphasK):
        flux = js[:,i]
        fit_obj = fit_objects[i]
        pars = parameters[i]
        psd = df_psd.iloc[-1].loc[str(Ks[i])]
        Es = Esmu[i]
        ### Plot energy flux j(E, alphaK) for target alphaK
        plots_fp.energy_y_data('flux', energy, flux, alphaK, '', True, 1, 0)

        plots_fp.check_fit_energy_psd(fit_obj, pars, psd, energy, flux,
                 fit_info, energy_range, Es, func, alphaK, unit_c, unit_flux)


def psd_lstar(df_psd, df_lstar, K, mu):
    lstar = df_lstar[str(K)][str(mu)].tolist()
    y = df_psd[str(K)][str(mu)].tolist()
    psd = u.Quantity(y)

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)

    ax.plot(lstar, psd, '-o', markersize= 4, color = 'b')

    title = 'K = ' + str(K) + ',    $\\mu$ = ' + str(mu)
    ax.set_yscale('log')
    ax.set_ylabel(f'PSD [{psd.unit}]', fontsize=17)
    ax.set_xlabel('$L^*$', fontsize=17)
#    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
    ax.tick_params(axis='both', labelsize=16)
    ax.set_title(title, fontsize = 20)

    plt.show()


def psd_time(df_psd, K, mu):
    time = df_psd[str(K)][str(mu)].index
    y = df_psd[str(K)][str(mu)].tolist()
    psd = u.Quantity(y)

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)

    ax.plot(time, psd, '-o', markersize= 4, color = 'b')

    title = 'K = ' + str(K) + ',    $\\mu$ = ' + str(mu)
    ax.set_yscale('log')
    ax.set_ylabel(f'PSD [{psd.unit}]', fontsize=17)
    ax.set_xlabel('Time', fontsize=17)
#    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
    ax.tick_params(axis='both', labelsize=16)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval = 1))  # Localizador de días
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))  # Formato de fecha
    ax.set_title(title, fontsize = 20)

    plt.show()
