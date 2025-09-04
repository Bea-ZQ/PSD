#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from astropy import units as u
from astropy import constants as const
import scipy.interpolate as scipy_int


'''
###############################################################################
### En este script están todas las funciones para hacer los cálculos de
### invariantes adiabáticas, necesarios para el cálculo de PSD
################################################################################
'''

######################### Segunda invariante K #################################


def K_mirror_points(alphas, model, inputs, unit_K):
    # Opción 1: Cálculo de K(alpha) usando find_mirror_points + make_lstar
    x_input, mag_input = inputs
    Ks  = []
    bs_mirr = []
    alphas =alphas.value

    # Para cada alfa, calculamos la posición del mirror point
    for i in range(len(alphas)):
    #    print('Alpha: ', alphas[i])
        dict_mirror = model.find_mirror_point(x_input, mag_input, alphas[i])
        posit_mirror = dict_mirror['POSIT']
        bmin = dict_mirror['bmin']      # Campo magnético en el punto mirror
        bs_mirr.append(bmin)

        x_geo = posit_mirror[0]
        y_geo = posit_mirror[1]
        z_geo = posit_mirror[2]

        # Calculamos I para cada alfa usando la posición de su mirror point
        x_input_alpha = {'x1':x_geo, 'x2':y_geo, 'x3':z_geo, 'dateTime':x_input['dateTime']}
        dict_lstar = model.make_lstar(x_input_alpha, mag_input)
        I_alpha = dict_lstar['xj'][0]
        K_alpha = I_alpha*np.sqrt(bmin)   # Unidades de Re*(nT**(1/2))
        Ks.append(K_alpha)

    return Ks*unit_K, bs_mirr*u.nT


def K_drift_bounce(alphas, model, inputs, unit_K):
    '''
    Esta es la única fn que modifiqué para el caso con b<0
    '''
    # Opción 2. Cálculo de K(alpha) usando drift_bounce_orbit
    x_input, mag_input = inputs
    Ks  = []
    bs_mirr = []
    bs_mirr2 = []
    alphas =alphas.value

    R0=1
    for i in range(len(alphas)):
#        print('Alpha: ', alphas[i])
        dict_bounce = model.drift_bounce_orbit(x_input, mag_input, alphas[i], R0)
        bmirr = dict_bounce['bmirr']
        bs_mirr.append(bmirr)

#        dict_mirror = model.find_mirror_point(x_input, mag_input, alphas[i])
#        bmirr2 = dict_mirror['bmin']      # Campo magnético en el punto mirror
#        bs_mirr2.append(bmirr2)

        #print('b_mirr: ', bmirr)
        I_alpha = dict_bounce['xj']
        #print('I: ', I_alpha)
        if bmirr >=0:
            K_alpha = I_alpha*np.sqrt(bmirr)  # Unidades de Re*(nT**(1/2))
        else:
            print('bmirr is negativo')
            K_alpha = np.nan
            print(K_alpha)
        #print('K: ', K_alpha)
        Ks.append(K_alpha)
        #print('----------------')
    return Ks*unit_K, bs_mirr*u.nT


def K_integration(alphas, model, inputs, unit_K):
    # Opción 3. Cálculo de K(alpha) usando mi propia integración numérica
    # y la información de la línea de campo
    x_input, mag_input = inputs
    Ks  = []
    bs_mirr = []
    R0 = 1
    alphas =alphas.value
    for i in range(len(alphas)):
        dict_fl = model.trace_field_line(x_input, mag_input, R0)
        posit_fl = dict_fl['POSIT']
        blocal = dict_fl['blocal']

        dict_mirror = model.find_mirror_point(x_input, mag_input, alphas[i])
        bmin = dict_mirror['bmin']
        bs_mirr.append(bmin)

        index_min = np.argmin(blocal)
        b1 = blocal[:index_min]
        b2 = blocal[index_min:]

        i1 = index_min - np.searchsorted(b1[::-1], bmin)
        i2 = index_min + np.searchsorted(b2, bmin)

        int_b_mag = blocal[i1:i2]
        int_posit_fl = posit_fl[i1:i2]

        x_fl = int_posit_fl[:,0]
        y_fl = int_posit_fl[:,1]
        z_fl = int_posit_fl[:,2]
        ds = np.sqrt(np.diff(x_fl)**2 + np.diff(y_fl)**2 + np.diff(z_fl)**2)

        ks = np.sqrt(bmin - int_b_mag)
        ks_avg = (ks[:-1] + ks[1:])*0.5

        line_integral = np.sum(ks_avg * ds) # Unidades de Re*(nT**(1/2))

        Ks.append(line_integral)

    return Ks*unit_K, bs_mirr*u.nT


def info_calculate_K(flag_K):
    if flag_K == '1':
        func = K_mirror_points
        str = 'find_mirror_points + make_lstar'
    elif flag_K == '2':
        func = K_drift_bounce
        str = 'drift_bounce_orbit'
    elif flag_K == '3':
        func = K_integration
        str = 'Integration'
    return func, str, flag_K


def interpolator_alpha(alphas, Ks):
    try:
        # k=1 especifica una interpolación lineal
        spline = scipy_int.make_interp_spline(np.log10(Ks.value[::-1]), alphas.value[::-1], k=1)
    except IndexError:
        print('spline not working')
        spline = False

    return spline


def interpolate_alpha_K(K, K_max, spline, flag=False):
    try:
        alpha = spline(np.log10(K.value))
    except TypeError:
        alpha = np.nan

    if (alpha<2) | (alpha >90):
        alpha = np.nan

    if K > K_max:
        alpha = np.nan

    if flag:
        return alpha*u.degree, (180.-alpha)*u.degree
    else:
        return alpha*u.degree


def calculate_mu(E, B, alpha):
    c2m0 = (const.m_e*const.c**2).to('MeV')
    F0 = 1 + (0.5*E)/c2m0
    mu = ((E*np.sin(alpha)**2)/B)*F0
    return mu


def calculate_E(mu, B, alpha):
    # B: magnitude magnetic field
    c2m0 = (const.m_e*const.c**2).to('MeV')
    F0 = 1 + ((2*B*mu)/(c2m0*(np.sin(alpha)**2)))
    E = c2m0*(-1 + np.sqrt(F0))
    return E


def Lstar_mirror_points(alpha, model, inputs):
    # Opción 1: Cálculo de Lstar) usando find_mirror_points + make_lstar
    # make l star calcula esta valor para partículas que tienen su punto de mirror
    # en la input position (90° PA). Yo quiero L star para partículas que en
    # la input position tengan PA alpha. Entonces lo que hago es encontrar la posición
    # en la cual van a rebotar las partículas con PA alpha (en tal posición van a
    # tener PA 90!) y usar esa posición para make_lstar

    ### Calculamos la posición del mirror point para partículas con PA alpha
    x_input, mag_input = inputs
    dict_mirror = model.find_mirror_point(x_input, mag_input, alpha.value)
    posit_mirror = dict_mirror['POSIT']

    x_geo = posit_mirror[0]
    y_geo = posit_mirror[1]
    z_geo = posit_mirror[2]

    ### Calculamos Lstar para partículas con PA apha_K usando la posición de su mirror point
    x_input_alpha = {'x1':x_geo, 'x2':y_geo, 'x3':z_geo, 'dateTime':x_input['dateTime']}
    dict_lstar = model.make_lstar(x_input_alpha, mag_input)
    lstar =  dict_lstar['Lstar'][0]
    return lstar


def Lstar_drift_bounce(alpha, model, inputs):
    # Opción 2. Cálculo de Lstar usando drift_bounce_orbit, que permite usar
    # un PA distinto de 90.
    x_input, mag_input = inputs
    R0=1
    dict_bounce = model.drift_bounce_orbit(x_input, mag_input, alpha.value, R0)
    lstar = dict_bounce['lstar']
    return lstar


def info_calculate_Lstar(flag_Lstar):
    if flag_Lstar == '1':
        func = Lstar_mirror_points
        str = 'find_mirror_points + make_lstar'
    elif flag_Lstar == '2':
        func = Lstar_drift_bounce
        str = 'drift_bounce_orbit'
    return func, str, flag_Lstar


def get_B_model(model, inputs):
    ### Obtenemos el campo magnético local en la posición del spacecraft
    ### usando IRBEM, recordar que está en nanoteslas
    print('* Using model for local magnenitc field calculation:')

    dict_mag_field = model.get_field_multi(*inputs)
    b_mag = dict_mag_field['Bl'][0]*u.nT
    return b_mag
