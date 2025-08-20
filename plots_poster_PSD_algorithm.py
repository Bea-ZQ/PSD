#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
import datetime
import time
from astropy import constants as const
from scipy.signal import argrelextrema

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u

def step2(list_Ks, list_alphas, list_splines, K, alpha, K_min, K_max, K_lim=0, save=0):

    fig, ax = plt.subplots()
    fig.set_size_inches(7, 6)
    label1 = 'IRBEM calculated $(K_i, \\alpha_i)$'
    label2 = 'Linear interpolation'

    if len(list_splines) == 1:
        Ks = list_Ks[0]
        alphas = list_alphas[0]
        print(Ks)
        print(alphas)
        ax.plot(Ks.value, alphas.value, 'o', color = 'b', label = label1)

        spline = list_splines[0]
        Ks_for_int = np.geomspace(K_min, K_max, 1000)
        alphas_interp = spline(np.log10(Ks_for_int.value))
        ax.plot(Ks_for_int.value, alphas_interp, color = 'r', label = label2)

    else:
        Ks_lin = list_Ks[0]
        alphas_lin = list_alphas[0]
        ax.plot(Ks_lin.value, alphas_lin.value, 'o', color = 'b', label = label1)

        Ks_geom = list_Ks[1]
        alphas_geom = list_alphas[1]
        ax.plot(Ks_geom.value, alphas_geom.value, 'o', color = 'b')

        spline_lin, spline_geom = list_splines
        Ks_lin_for_int = np.linspace(K_min, K_lim, 1000)
        Ks_geom_for_int = np.linspace(K_lim, K_max, 1000)
        alphas_lin_interp = spline_lin(np.log10(Ks_lin_for_int.value))
        alphas_geom_interp = spline_geom(np.log10(Ks_geom_for_int.value))

        ax.plot(Ks_lin_for_int.value, alphas_lin_interp, color = 'r', label = label2)
        ax.plot(Ks_geom_for_int.value, alphas_geom_interp, color='r')

    ax.set_xlim([8,3000])
    ax.set_ylim([0,95])
    ax.set_xlabel('K [%s]' % str(K.unit), fontsize =17)
    ax.set_ylabel('Pitch angle $\\alpha$ [°]',fontsize = 17)
    ax.set_xscale('log')
    ax.tick_params(axis='both', labelsize=16)
    ax.hlines(y=alpha.value, xmin = 8, xmax = K.value, color= 'r')
    ax.vlines(x=K.value,ymin =0,  ymax = alpha.value, color= 'r')

    ax.plot(K.value, alpha.value, 'o', color = 'r', label = 'Target $(K,\\alpha_K)$')
    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=14)
    if save == 1:
        fig.savefig('plot1.pdf')
    plt.show()


def step5(fit_info, fedu, list_bins, fit_results, js_alphaK, alphaK, i, save=0):

    fit_opt, title_plot, func, _, _ = fit_info
    energy_bins, alpha_bins = list_bins
    fit_flag, parameters, max_values = fit_results

    alphas = np.linspace(20,160, 1000)*u.degree
    alphas_rad = np.radians(alphas)


    colors = ['b', 'g']
    fig, ax = plt.subplots()
    fig.set_size_inches(7, 6)
#    ax.set_aspect('equal')
    if fit_flag[i]:
#        title_plot = 'Fitted data. Function ' + func_opt
        par_opt = parameters[i]
        max_val = max_values[i]
        j = func(alphas_rad, *par_opt)
        ax.plot(alphas, j*max_val, color = 'r', label = 'Fitted curve')
        ax.set_title(title_plot)
        ax.hlines(y=js_alphaK[i].value, xmin = 20, xmax = alphaK.value, color= 'r')
        ax.vlines(x=alphaK.value, ymin =j[0]*max_val.value,  ymax = js_alphaK[i].value, color= 'r')
        ax.set_ylim([j[0]*max_val.value, None])

    else:
        ax.set_title('No fit')

    ax.set_ylabel('Electron flux [%s]' % str(fedu.unit), fontsize =17)
    ax.set_xlabel('Pitch angle $\\alpha$ [°]', fontsize =17)

    flux = fedu[:,i]
    ax.plot(alpha_bins, flux, 'o', color = 'b', markersize = 10, label = 'Flux data for %s' %str(energy_bins[i]))
    ax.plot(alphaK, js_alphaK[i], 'o', color = 'r', label = 'Flux at target $\\alpha_K$')

    ax.set_title(title_plot, fontsize = 20)
    ax.set_xlim([20,160])
    ax.tick_params(axis='both', labelsize=16)

    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=14)
    ax.set_yscale('log')
    if save ==1:
        fig.savefig('plot2.pdf')
    plt.show()


def step6(fit_results, energy_bins, y, fit_info, energy_range, Emu, save=0):
    func, parms, y_Emu = fit_results
    fit_opt, str_func, _,  _, _ = fit_info
    E_min, E_max = energy_range

    energies = np.linspace(E_min,E_max, 1000).to(u.MeV)
    log_energies = np.log10(energies.value)

    fig, ax = plt.subplots()
    fig.set_size_inches(7, 6)
    ax.set_yscale('log')

    fit = isinstance(func, bool)
    print(fit)
    if fit:
        print('in if')
        return

    if fit_opt[1] == 'flux':
        print('in if flux')
        if fit_opt[0] == 'pl':
            label_fit = 'Power law fit'
            #'$j(E) = aE^{-k}$'
        elif fit_opt[0] == 'spl':
            label_fit = 'Spline'
        elif fit_opt[0] == 'lin':
            label_fit = 'Linear interpolation '

        log_j = func(log_energies, *parms)
        j = 10**log_j*y.unit
        ax.plot(energies, j, color = 'r', label = label_fit)
        ax.set_ylabel('Electron flux [%s]' % str(y.unit), fontsize =17)
        label1 = 'Flux at target $\\alpha_K$'
        label2 = 'Flux at target $\\alpha_K$ and $E_\\mu$'

    elif fit_opt[1] == 'psd':
        if fit_opt[0] == 'exp':
            label_fit = 'Exponential fit'
            f = func(energies.value, *parms)*y.unit
            ax.plot(energies, f, color = 'r', label = label_fit)
        elif fit_opt[0] == 'spl':
            log_f = func(log_energies, *parms)
            f = 10**log_f*y.unit
            ax.plot(energies, f, color = 'r', label = 'Spline')
        elif fit_opt[0] == 'lin':
            log_f = func(log_energies, *parms)
            f = 10**log_f*y.unit
            ax.plot(energies, f, color = 'r', label = 'Linear interpolation')
        ax.set_ylabel('Electron PSD [%s]' % str(y.unit), fontsize =17)
        label1 = 'PSD at target $\\alpha_K$'
        label2 = 'PSD at target $E_\\mu'

    ax.plot(energy_bins,y.value, '*', color = 'b', markersize = 10, label = label1)
    ax.plot(Emu, y_Emu, 'o', color = 'r', label = label2)
    ax.set_xlabel('Energy [%s]' % str(energy_bins.unit), fontsize=17)
    ax.tick_params(axis='both', labelsize=16)
    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=14)
    ax.set_title(str_func, fontsize = 20)
    ax.set_xscale('log')
    ax.set_xlim([E_min.value, E_max.value])
    if save ==1:
        fig.savefig('plot3.pdf')
    plt.show()
    plt.close(fig)
