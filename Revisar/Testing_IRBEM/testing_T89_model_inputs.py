###############################################################################
### En este script estoy creando usando magnétic inputs y x inputs que provienen
### de los datos REPT y OMNI para el modelo T89
################################################################################

import os, sys
import datetime
import pandas as pd
import numpy as np
from inputs_utils import get_GEO_coordinates, dict_x_input
from inputs_utils import T89_magnetic_inputs, filter_dates, dict_magnetic_input_T89
from plot_utils import plot_orbit, plot_field_line, plot_drift_shell
import matplotlib.pyplot as plt
import IRBEM

################################################################################
### Definiendo directorios importantes
################################################################################


home_path = os.path.expanduser("~")
rb_path = '/Work/Research/Radiation_Belts'
data_rept_path = '/data/rbsp/'
data_omni_path = '/data/omni/'
codes_path = '/RB_codes/Data'

sys.path.append(home_path + rb_path + codes_path)
from rept_data_utils import load_CDFfiles_REPT, download_CDFfiles_REPT
from omni_data_utils import load_CDFfiles_OMNI, download_CDFfiles_OMNI
from omni_data_utils import get_local_dir_OMNI, get_remote_dir_OMNI, get_filename_OMNI
from omni_data_utils import read_CDFfile_OMNI

remote_omni_dir = "https://cdaweb.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/"
data_omni_dir = home_path + rb_path + data_omni_path

remote_rept_dir = "https://rbsp-ect.newmexicoconsortium.org/data_pub/"
data_rept_dir = home_path + rb_path + data_rept_path



################################################################################
###         Definiendo fechas REPT y OMNI para el análisis: CREAR UNA FUNCIÓN QUE HAGA ESTO BIEN
################################################################################
sdate = '2013-05-30'
edate = '2013-05-31'

flag = 7 > int(sdate[5:7])
if flag:
    sdate_omni = sdate[:5] + '01-01'
else:
    sdate_omni = sdate[:5] + '07-01'
edate_omni = pd.to_datetime(edate) + pd.Timedelta(days=1)



################################################################################
###                     Descarga de datos rept y omni                        ###
################################################################################

#download_CDFfiles_REPT(sdate, edate, data_rept_dir, remote_rept_dir, probe='both', level="3")
#download_CDFfiles_OMNI(sdate_omni, edate_omni, remote_omni_dir, data_omni_dir, '1h', '')


################################################################################
###                       Lectura de datos rept y omni                       ###
################################################################################

# LAS POSICIONES DE REPT ESTÁN EN GEO, EN KM, TENGO QUE NORMALIZAR POR EL RADIO DE LA TIERRA
relevant_keys_rept = ['Epoch', 'Position']
Re = 6371

[probeA, probeB]= load_CDFfiles_REPT(sdate, edate, data_rept_dir, relevant_keys_rept, 'both', '3', 'fedu')
[reptA_info, feduA, reptA_metadata, feduA_metadata] = probeA
[reptB_info, feduB, reptB_metadata, feduB_metadata] = probeB

omni1h_info, omni1h_metadata = load_CDFfiles_OMNI(sdate_omni, edate_omni, data_omni_dir, '1h', '')


################################################################################
###              Recuperamos tiempo y posiciciones de datos REPT             ###
################################################################################

posit_A = reptA_info[['Epoch', 'posit_x1', 'posit_x2', 'posit_x3']]
posit_B = reptB_info[['Epoch', 'posit_x1', 'posit_x2', 'posit_x3']]

dates_A, x_geo_A, y_geo_A, z_geo_A = get_GEO_coordinates(posit_A, Re)
dates_B, x_geo_B, y_geo_B, z_geo_B = get_GEO_coordinates(posit_B, Re)


################################################################################
#                        Creamos diccionario x_input
################################################################################

x_input_A = dict_x_input(dates_A, x_geo_A, y_geo_A, z_geo_A)
x_input_B = dict_x_input(dates_B, x_geo_B, y_geo_B, z_geo_B)

x_input_A_float = dict_x_input(dates_A[0], x_geo_A[0], y_geo_A[0], z_geo_A[0])
x_input_B_float = dict_x_input(dates_B[0], x_geo_B[0], y_geo_B[0], z_geo_B[0])

################################################################################
################################################################################

df_omni, metadata_omniT89 = T89_magnetic_inputs(omni1h_info, omni1h_metadata)
df_omni_filt = filter_dates(df_omni, sdate, edate)


################################################################################
#                  Creamos diccionario magnetic input for T89
################################################################################

#.loc es el índice
#.iloc es la posición



mag_input_A = dict_magnetic_input_T89(dates_A, df_omni_filt)
mag_input_B = dict_magnetic_input_T89(dates_B, df_omni_filt)

mag_input_A_float = dict_magnetic_input_T89(dates_A.iloc[4367], df_omni_filt, True)
mag_input_B_float = dict_magnetic_input_T89(dates_B.iloc[12854], df_omni_filt, True)

################################################################################
# Probando el modelo: trace field line y get_field_multi
################################################################################

options = [0,0,0,0,0]
kext_t89 = 'T89'
model_t89 = IRBEM.MagFields(options=options, kext=kext_t89, verbose = True, sysaxes=1)


# get field multi
out_dict = model_t89.get_field_multi(x_input_B, mag_input_B)

bx_geo = out_dict['BxGEO']
by_geo = out_dict['ByGEO']
bz_geo = out_dict['BzGEO']
b_mag = out_dict['Bl']

# Trace field line
R0 = 1
dens = 2
out_dict_float = model_t89.trace_field_line(x_input_B_float, mag_input_B_float, R0)
posit = out_dict_float['POSIT']

# Now plot the field lines
ax = plot_field_line(posit, dens)
ax.scatter(x_input_B_float['x1'], x_input_B_float['x2'], x_input_B_float['x3'], color='red', s=100, label='initial pos')
ax.legend()
plt.show()

################################################################################
# Probando el modelo: drift bounce orbit
################################################################################

alpha = 60
R0 = 1
out_dict_float = model_t89.drift_bounce_orbit(x_input_B_float, mag_input_B_float, alpha, R0)

posit = out_dict_float['POSIT']
Nposit = out_dict_float['Nposit']

dens = 5
ax = plot_drift_shell(posit, Nposit, dens, 2)
ax.scatter(x_input_B_float['x1'], x_input_B_float['x2'], x_input_B_float['x3'], color='red', s=100, label='initial pos')
ax.legend()
plt.show()
