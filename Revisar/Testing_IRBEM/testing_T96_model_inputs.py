###############################################################################
### En este script estoy creando usando magnétic inputs y x inputs que provienen
### de los datos REPT y OMNI para el modelo T96
################################################################################

import os, sys
import datetime
import pandas as pd
import numpy as np
from plot_inputs import plot_orbit, plot_field_line, plot_drift_shell
from inputs_utils import get_GEO_coordinates, dict_x_input
from inputs_utils import T96_magnetic_inputs, filter_dates, dict_mag_input_T96, check
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
sdate = '2013-05-31'
edate = '2013-05-31'

sdate_omni5m = sdate[:8] + '01'
edate_omni5m = pd.to_datetime(edate) + pd.Timedelta(days=1)

flag = 7 > int(sdate[5:7])
if flag:
    sdate_omni1h = sdate[:5] + '01-01'
else:
    sdate_omni1h = sdate[:5] + '07-01'
edate_omni1h = pd.to_datetime(edate) + pd.Timedelta(days=1)


file_type = 'hro2'

################################################################################
###                     Descarga de datos rept y omni                        ###
################################################################################

#download_CDFfiles_REPT(sdate, edate, data_rept_dir, remote_rept_dir, probe='both', level="3")
#download_CDFfiles_OMNI(sdate_omni, edate_omni, remote_omni_dir, data_omni_dir, '5min', file_type)

################################################################################
###                       Lectura de datos rept y omni                       ###
################################################################################

# LAS POSICIONES DE REPT ESTÁN EN GEO, EN KM, TENGO QUE NORMALIZAR POR EL RADIO DE LA TIERRA
relevant_keys_rept = ['Epoch','MLAT', 'MLT', 'B_calc', 'I', 'L', 'L_star', 'Position']
Re = 6371

[probeA, probeB]= load_CDFfiles_REPT(sdate, edate, data_rept_dir, relevant_keys_rept, 'both', '3', 'fedu')
[reptA_info, feduA, reptA_metadata, feduA_metadata] = probeA
[reptB_info, feduB, reptB_metadata, feduB_metadata] = probeB

omni5m_info, omni5m_metadata = load_CDFfiles_OMNI(sdate_omni5m, edate_omni5m, data_omni_dir, '5min', file_type)
omni1h_info, omni1h_metadata = load_CDFfiles_OMNI(sdate_omni1h, edate_omni1h, data_omni_dir, '1h', '')


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


df_omni1h, metadata_omni1h = T96_magnetic_inputs([omni1h_info, omni1h_metadata], '', '1h')
df_omni5m, metadata_omni5m = T96_magnetic_inputs([omni1h_info, omni1h_metadata],
                                                [omni5m_info, omni5m_metadata], '5min')

df_omni1h_filt = filter_dates(df_omni1h, sdate, edate)
df_omni5m_filt = filter_dates(df_omni5m, sdate, edate)



################################################################################
#                  Creamos diccionario magnetic input for T96
################################################################################

#.loc es el índice
#.iloc es la posición

mag_input1h_A = dict_mag_input_T96(dates_A, df_omni1h_filt)
mag_input1h_B = dict_mag_input_T96(dates_B, df_omni1h_filt)

# Checked
mag_input1h_A_float = dict_mag_input_T96(dates_A.iloc[3367], df_omni1h_filt, True)
mag_input1h_B_float = dict_mag_input_T96(dates_B.iloc[7232], df_omni1h_filt, True)

mag_input5m_A = dict_mag_input_T96(dates_A, df_omni5m_filt)
mag_input5m_B = dict_mag_input_T96(dates_B, df_omni5m_filt)

# Checked
mag_input5m_A_float = dict_mag_input_T96(dates_A.iloc[3367], df_omni5m_filt, True)
mag_input5m_B_float = dict_mag_input_T96(dates_B.iloc[7232], df_omni5m_filt, True)


################################################################################
# Probando modelo T96 array
################################################################################
options = [0,0,0,0,0]
kext = 'T96'
model = IRBEM.MagFields(options=options, kext=kext, verbose = True, sysaxes=1)


out_dict1h = model.get_field_multi(x_input_B, mag_input1h_B)

#bx_geo = out_dict1h['BxGEO']
#by_geo = out_dict1h['ByGEO']
#bz_geo = out_dict1h['BzGEO']
#b_mag = out_dict1h['Bl']

out_dict5m = model.get_field_multi(x_input_A, mag_input5m_A)





'''
################################################################################
# Probando modelo T96 float
################################################################################

R0 = 1
dens = 2
#out_dict1h_float = model.trace_field_line(x_input_B_float, mag_input1h_B_float, R0)
#posit1h = out_dict1h_float['POSIT']

# Now plot the field lines
#ax = plot_field_line(posit1h, dens)
#ax.scatter(x_input_B_float['x1'], x_input_B_float['x2'], x_input_B_float['x3'], color='red', s=100, label='initial pos')
#ax.legend()
#plt.show()


out_dict5m_float = model.trace_field_line(x_input_A_float, mag_input5m_A_float, R0)
posit5m = out_dict5m_float['POSIT']

# Now plot the field lines
ax = plot_field_line(posit5m, dens)
ax.scatter(x_input_A_float['x1'], x_input_A_float['x2'], x_input_A_float['x3'], color='red', s=100, label='initial pos')
ax.legend()
plt.show()
'''

check(dates_A, mag_input5m_A)
print(df_omni5m_filt.iloc[-50:])
