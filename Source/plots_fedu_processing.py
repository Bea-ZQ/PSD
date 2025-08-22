#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
from astropy import units as u
#matplotlib.rcParams['text.usetex'] = True
#plt.rcParams['font.family'] = 'serif'


'''
###############################################################################
### En este script están todas las funciones para visualizar los ajustes
### de los flujos de electrones, necesarios para el cálculo de PSD
################################################################################
'''

#################################### PA fit ####################################

def PA_flux_data(energy_bins, alpha_bins, fedu, instrument, save_path, prefix, show=1, save=1):
    for i in range(1, len(energy_bins)+1):
        save_name_i = os.path.join(save_path, 'PA_flux', '%s_%s_channel_%d.png' % (instrument,prefix, i))
        fig_i, ax_i = plt.subplots()
    #    print('----')
    #    print(i)
        flux = fedu[:,i-1]
#        print(len(flux))
        ax_i.plot(alpha_bins, flux, '*', label = energy_bins[i-1], color = 'b')

#        ax_i.plot(alpha_bins, flux, '*',  color = 'b')

        ax_i.set_ylabel('Flux [%s]' % str(fedu.unit))
        ax_i.set_xlabel('Pitch angle [°]')
        ax_i.set_yscale('log')
        ax_i.legend()
        ax_i.set_title('%s electron fluxes\nChannel %d' %(instrument, i))
        if save ==1:
            fig_i.savefig(save_name_i)
        if show ==1:
            plt.show()
        plt.close(fig_i)
    return


def PA_flux_function(fit_opt, str_title, func, parms, max_value, save_path,
                     show=1, save=1):

    alphas = np.linspace(20,160, 1000)*u.degree
    alphas_rad = np.radians(alphas)

    j = func(alphas_rad, *parms)

    title = 'Function %s to fit j(alpha)\n%s' % (fit_opt, str_title)

    fig, ax = plt.subplots()
    ax.plot(alphas, j*max_value, color = 'r')
    ax.set_ylabel('Flux [%s]' % str(max_value.unit))
    ax.set_xlabel('Pitch angle [°]')
    ax.set_yscale('log')
    ax.set_title(title)
    if save ==1:
        save_name = os.path.join(save_path, 'PA_flux', 'Function%s_to_fit_PA_flux.png' % fit_opt)
        fig.savefig(save_name)
    if show ==1:
        plt.show()
    return fig, ax


def check_fit_PA_flux(fit_opts, list_func, fedu, list_bins, fit_results, js_alphaK, alphaK,
    instrument, save_path, show=1, save=1):

    energy_bins, alpha_bins = list_bins
    fit_flag, parameters, max_values = fit_results[:3]
    N= len(energy_bins)
    j=0
    for i in range(N):
        if fit_flag[i]:
            title = '%s channel %d fitted data' %(instrument, i+1)
            par_opt = parameters[i]
            max_val = max_values[i]
            opt = fit_opts[i]
            func = list_func[i]
            fig, ax = PA_flux_function(opt, title, func, par_opt, max_val, '', 0, 0)
            ax.plot(alphaK, js_alphaK[j], 'o', color = 'r', label = 'Fitted flux')
            j+=1
        else:
            opt = fit_opts[i]
            fig, ax = plt.subplots()
            title = 'Function %s to fit j(alpha)\n%s channel %d: No fit' % (opt, instrument, i+1)
            ax.set_title(title)
            ax.set_ylabel('Flux [%s]' % str(fedu.unit))
            ax.set_xlabel('Pitch angle [%s]' % str(alpha_bins.unit))

        flux = fedu[:,i]
        ax.plot(alpha_bins, flux, '*', color = 'b', label = energy_bins[i])
        ax.legend()
        ax.set_yscale('log')
        if save ==1:
            save_name = os.path.join(save_path, 'PA_flux', 'Fit_%s_function%s_channel%d.png' % (instrument, opt, i+1))
            fig.savefig(save_name)
        if show ==1:
            plt.show()
        plt.close(fig)
    return



################################## Energy fit ##################################

def energy_y_data(data_type, energy, y_data, alphaK, save_path, flag_log=True, show=1, save=1):
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 6)
    ax.set_xlabel('Energy [%s]' % str(energy.unit))
    ax.set_yscale('log')

    if data_type == 'flux':
        label = 'Interpolated j(E, alphaK)'
        ax.set_ylabel('Flux [%s]' % str(y_data.unit))
        ax.set_title('Electron fluxes at alphaK = %0.2f°' %alphaK.value)
    elif data_type == 'psd':
        label = 'Interpolated psd(E, alphaK)'
        ax.set_ylabel('PSD [%s]' % str(y_data.unit))
        ax.set_title('Electron PSD at alphaK = %0.2f°' %alphaK.value)

    ax.plot(energy, y_data, '*', color = 'r', label = label)

    if flag_log:
        ax.set_xscale('log')
    ax.legend()
    if save ==1:
        save_name = os.path.join(save_path, 'Energy_data', 'Energy_%s_alphaK%.2f.png' %(data_type, alphaK.value))
        fig.savefig(save_name)
    if show ==1:
        plt.show()
    plt.close(fig)
    return


def check_fit_energy_data(fit_results, energy_bins, y, fit_info,
    energy_range, Emu, save_path, show=1, save=1):
    func, parms, y_Emu = fit_results
    fit_opt, str_func, _, _, _ = fit_info
    E_min, E_max = energy_range

    energies = np.linspace(E_min,E_max, 1000).to(u.MeV)
    log_energies = np.log10(energies.value)

    fit = not isinstance(func, bool)

    if fit:
        label_title = ''
    else:
        label_title = ': No fit'

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 6)

    if fit_opt[1] == 'flux':
        print('\nChecking flux fit')
        label = 'Flux at target Emu'
        if fit:
            log_j = func(log_energies, *parms)
            j = 10**log_j*y.unit
            ax.plot(energies, j, color = 'r')
        ax.set_ylabel('Flux [%s]' % str(y.unit))
    elif fit_opt[1] == 'psd':
        print('\nPSD fitting')
        label = 'PSD at target Emu'

        if fit_opt[0] == 'exp':
            if fit:
                f = func(energies.value, *parms)*y.unit
                ax.plot(energies, f, color = 'r')
        elif (fit_opt[0] == 'spl') | (fit_opt[0] == 'lin'):
            if fit:
                log_f = func(log_energies, *parms)
                f = 10**log_f*y.unit
                ax.plot(energies, f, color = 'r')
        ax.set_ylabel('PSD [%s]' % str(y.unit))

    title = 'Function %s to fit %s\n%s%s' % (fit_opt[0], fit_opt[1], str_func, label_title)
    ax.plot(Emu, y_Emu, 'o', color = 'r', label = label)
    ax.plot(energy_bins.value, y.value, '*', color = 'b', label = 'Data to fit')

    ax.set_xlabel('Energy [%s]' % str(energies.unit))
    ax.set_title(title)
    if fit:
        ax.set_xscale('log')
        ax.set_yscale('log')
    if save ==1:
        save_name = os.path.join(save_path, 'Energy_data', 'Fit_%s_function_%s_Emu%.2f.png' % (fit_opt[1], fit_opt[0], Emu.value))
        fig.savefig(save_name)
    if show ==1:
        plt.show()
    plt.close(fig)
    return


def check_fit_energy_psd(func, parms, f_Emu, energy_bins, flux_data, fit_info,
    energy_range, Emu, psd_func, alphaK, units_psd, show=1):
    fit_opt, str_func, _, _, _ = fit_info
    E_min, E_max = energy_range

    energies = np.linspace(E_min,E_max, 1000).to(u.MeV)
    log_energies = np.log10(energies.value)

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 6)
    label = 'PSD at target Emu'
    psd_data = psd_func(energy_bins, flux_data, units_psd)
    fit = isinstance(func, bool)
    print(fit)
    if fit:
        print('in if')
        label_title = ': No fit'
    else:
        label_title = ''
        if fit_opt[1] == 'flux':
            print('\nChecking flux fit')
            log_j = func(log_energies, *parms)
            j = 10**log_j*flux_data.unit
            f = psd_func(energies, j, units_psd)
        elif fit_opt[1] == 'psd':
            print('\nPSD fitting')
            if fit_opt[0] == 'exp':
                f = func(energies.value, *parms)*psd_data.unit
            elif (fit_opt[0] == 'spl') | (fit_opt[0] == 'lin'):
                log_f = func(log_energies, *parms)
                f = 10**log_f*psd_data.unit
        ax.plot(energies, f, color = 'r')
        ax.plot(energy_bins.value, psd_data.value, '*', color = 'b', label = 'Data to fit')
        for j in range(len(Emu)):
            if j==0:
                ax.plot(Emu[j], f_Emu.iloc[j], 'o', color = 'r', label = label)
            else:
                ax.plot(Emu[j], f_Emu.iloc[j], 'o', color = 'r')

    title = 'Function %s to fit %s\n%s%s\nalphaK = %0.2f°' % (fit_opt[0], fit_opt[1], str_func, label_title, alphaK.value)
    ax.set_ylabel('PSD [%s]' % str(psd_data.unit))
    ax.set_xlabel('Energy [%s]' % str(energies.unit))
    ax.set_title(title)
    ax.set_xscale('log')
    ax.set_yscale('log')

    if show ==1:
        plt.show()
    plt.close(fig)
    return
