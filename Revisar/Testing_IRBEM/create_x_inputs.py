###############################################################################
### En este script estoy creando X inputs para los modelos, para ver si
### las funciones están haciendo un buen trabajo
################################################################################

import os, sys
import datetime
import pandas as pd
import numpy as np
from create_inputs_utils import get_GEO_coordinates, dict_x_input
from plot_utils import plot_orbit
import IRBEM
'''
################################################################################
### Definiendo directorios importantes
################################################################################
'''

home_path = os.path.expanduser("~")
rb_path = '/Work/Research/Radiation_Belts'
data_path = '/data/rbsp/'
codes_path = '/RB_codes/Data'

sys.path.append(home_path + rb_path + codes_path)
from rept_data_utils import load_CDFfiles_REPT, download_CDFfiles_REPT

remote_dir = "https://cdaweb.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/"
data_dir = home_path + rb_path + data_path


################################################################################
### Definiendo fechas para el análisis
################################################################################
start_date = '2013-05-31'
end_date = '2013-05-31'
Re = 6371

################################################################################
###                           Descarga de datos                              ###
################################################################################

#download_CDFfiles_REPT(start_date, end_date, data_dir, remote_dir, probe='both', level="3")

################################################################################
###                             Lectura de datos                             ###
################################################################################

# LAS POSICIONES DE REPT ESTÁN EN GEO, EN KM, TENGO QUE NORMALIZAR POR EL RADIO DE LA TIERRA
relevant_keys = ['Epoch','MLAT', 'MLT', 'B_calc', 'I', 'L', 'L_star', 'Position']

#[probeA] = load_CDFfiles_REPT(start_date, end_date, data_dir, variables, 'a', '3', 'fedu')
#[probeB] = load_CDFfiles_REPT(start_date, end_date, data_dir, variables, 'b', '3', 'fedu')
[probeA, probeB]= load_CDFfiles_REPT(start_date, end_date, data_dir, relevant_keys, 'both', '3', 'fedu')

[reptA_info, feduA, reptA_metadata, feduA_metadata] = probeA
[reptB_info, feduB, reptB_metadata, feduB_metadata] = probeB

posit_A = reptA_info[['Epoch', 'posit_x1', 'posit_x2', 'posit_x3']]
posit_B = reptB_info[['Epoch', 'posit_x1', 'posit_x2', 'posit_x3']]

################################################################################
# Graficando posiciones
################################################################################


dates_A, x_geo_A, y_geo_A, z_geo_A = get_GEO_coordinates(posit_A, Re)
dates_B, x_geo_B, y_geo_B, z_geo_B = get_GEO_coordinates(posit_B, Re)

geo_array_A = np.array([x_geo_A, y_geo_A, z_geo_A]).T
geo_array_B = np.array([x_geo_B, y_geo_B, z_geo_B]).T

sm_array_A = IRBEM.Coords().transform(time=dates_A, pos=geo_array_A, sysaxesIn='GEO', sysaxesOut='SM')
sm_array_B = IRBEM.Coords().transform(time=dates_B, pos=geo_array_B, sysaxesIn='GEO', sysaxesOut='SM')

ax_lim = 6
#plot_orbit(geo_array_A, geo_array_B, 'GEO', ax_lim)
#plot_orbit(sm_array_A, sm_array_B, 'SM', ax_lim)


################################################################################
# creando diccionario input de los modelos
################################################################################

x_input_A = dict_x_input(dates_A, x_geo_A, y_geo_A, z_geo_A)
x_input_B = dict_x_input(dates_B, x_geo_B, y_geo_B, z_geo_B)

x_input_A_float = dict_x_input(dates_A[0], x_geo_A[0], y_geo_A[0], z_geo_A[0])
x_input_B_float = dict_x_input(dates_B[0], x_geo_B[0], y_geo_B[0], z_geo_B[0])
