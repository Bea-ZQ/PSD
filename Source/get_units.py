#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from astropy import units as u
from astropy import constants as const


'''
###############################################################################
### En este script están todas las funciones para obtener unidades
################################################################################
'''

######################### Segunda invariante K #################################

def adiab_inv():
    Re_m = const.R_earth
    Re_km = Re_m.to(u.km).value
    unit_Re = u.def_unit('Re', const.R_earth.value*const.R_earth.unit)
    unit_K = unit_Re*np.sqrt(1*u.nT)
    unit_mu = u.MeV/u.nT
    return Re_km, unit_Re, unit_K, unit_mu


def flux():
    unit_flux_mev =(u.cm**(2)*u.s*u.sr*u.MeV)**(-1)
    unit_flux_kev =(u.cm**(2)*u.s*u.sr*u.keV)**(-1)
    return unit_flux_mev, unit_flux_kev


def psd():
    unit_c = u.def_unit('c', const.c)
    unit_psd = (unit_c**3)/((u.MeV**3)*(u.cm**3))
    return unit_c, unit_psd
