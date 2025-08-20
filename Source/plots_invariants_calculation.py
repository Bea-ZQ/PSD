#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
#matplotlib.rcParams['text.usetex'] = True
#plt.rcParams['font.family'] = 'serif'


'''
###############################################################################
### En este script están todas las funciones para visualizar que los cálculos
### que involucran invariantes adiabáticas sean correctos
################################################################################
'''

############################### Segunda invariante #############################

def plot_alpha_K(list_Ks, list_alphas, list_labels, list_colors):
    fig, ax = plt.subplots()
    for i in range(len(list_Ks)):
        Ks = list_Ks[i]
        alphas = list_alphas[i]
        label = list_labels[i]
        color = list_colors[i]
        ax.plot(np.log10(Ks.value), alphas.value, '*', color = color, label = label)

    ax.set_ylabel('alphas')
    ax.set_xlabel('log10(K)')
    return ax


def interpolation_alpha_K(list_Ks, list_alphas, list_splines, K_min, K_max, K_lim = 0):
    # Esta función es para verificar visualmente que la interpolación esté funcionando

    list_colors = ['silver', 'b', 'r']
    if len(list_splines) == 1:
        list_labels = ['Calculated with MF model', 'Used for spline']
        spline = list_splines[0]
    else:
        list_labels = ['Calculated with MF model', 'Used for spline lin', 'Used for spline geom']
        spline_lin, spline = list_splines

    ax = plot_alpha_K(list_Ks, list_alphas, list_labels, list_colors)

    N = 100
    M = 5
    K_range = np.linspace(K_min, K_max, N)
    list_idx = np.random.randint(0, len(K_range), M)
    for idx in list_idx:
        K = K_range[idx]

        if K > K_lim:
            alpha = spline(np.log10(K.value))
        else:
            alpha = spline_lin(np.log10(K.value))
        ax.plot(np.log10(K.value), alpha, 'o', color = 'g')
    ax.legend()
    plt.show()


def calculation_E_mu(list_E, list_mu, E, mu):
    fig, ax = plt.subplots()
    ax.plot(list_E.value, list_mu.value, '*', color = 'b')
    ax.set_xlabel('Energy [%s]' % str(list_E.unit))
    ax.set_ylabel('mu [%s]' % str(list_mu.unit))
    ax.plot(E.value, mu, 'o', color = 'r', label = 'target')
    plt.show()
