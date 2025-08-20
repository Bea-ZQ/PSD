###############################################################################
### En este script estoy creando magnétic inputs para los modelos, para ver si
### las funciones están haciendo un buen trabajo
################################################################################

import os, sys
import datetime
import pandas as pd
import numpy as np
from create_inputs_utils import T89_magnetic_inputs, T96_magnetic_inputs, filter_dates
from plot_utils import plot_omni

'''
################################################################################
### Definiendo directorios importantes
################################################################################
'''

home_path = os.path.expanduser("~")
rb_path = '/Work/Research/Radiation_Belts'
data_path = '/data/omni/'
codes_path = '/RB_codes/Data'

sys.path.append(home_path + rb_path + codes_path)
from omni_data_utils import load_CDFfiles_OMNI, download_CDFfiles_OMNI
from omni_data_utils import get_local_dir_OMNI, get_remote_dir_OMNI, get_filename_OMNI
from omni_data_utils import read_CDFfile_OMNI

remote_dir = "https://cdaweb.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/"
data_dir = home_path + rb_path + data_path


################################################################################
### Definiendo fechas para el análisis
################################################################################
start_year = '2015'
end_year = '2015'
file_type = 'hro'

start_date = start_year + '-01-01'
end_date = end_year + '-12-31'

print('start_date: ', start_date)
print('end_date: ', end_date)

################################################################################
###                           Descarga de datos                              ###
################################################################################

#download_CDFfiles_OMNI(start_date, end_date, remote_dir, data_dir, '1h', '')
#download_CDFfiles_OMNI(start_date, end_date, remote_dir, data_dir, '5min', file_type)

################################################################################
###                             Lectura de datos                             ###
################################################################################

df_omni1h, dict_metadata_omni1h = load_CDFfiles_OMNI(start_date, end_date, data_dir, '1h', '')
df_omni5m, dict_metadata_omni5m = load_CDFfiles_OMNI(start_date, end_date, data_dir, '5min', file_type)


################################################################################
###                             Magnetic inputs                              ###
################################################################################

df_inputsT89, metadataT89 = T89_magnetic_inputs(df_omni1h, dict_metadata_omni1h)

df_inputsT96_1h, metadataT96_1h = T96_magnetic_inputs([df_omni1h, dict_metadata_omni1h], '', '1h')
df_inputsT96_5m, metadataT96_5m = T96_magnetic_inputs([df_omni1h, dict_metadata_omni1h], [df_omni5m, dict_metadata_omni5m], '5min')


################################################################################
###                             Filtramos fechas                             ###
################################################################################

sdate= '2015-03-01'
edate= '2015-03-31'


df_T96_plot = filter_dates(df_inputsT96_1h, sdate, edate)


################################################################################
###                               Graficamos                                 ###
################################################################################

x = df_T96_plot['epoch']
y_dst = df_T96_plot['Dst']
y_Bz = df_T96_plot['Bz']
y_By = df_T96_plot['By']
y_pdyn = df_T96_plot['Pdyn']

plot(x, y_dst, 'Dst')
plot(x, y_Bz, 'Bz')
plot(x, y_By, 'By')
plot(x, y_pdyn, 'Pdyn')




'''
date_array = pd.date_range(start=start_date, end=end_date, freq='6MS')
print('date_array', date_array)

filename = get_filename_OMNI(date_array[3], resolution, file_type)
local_dir = get_local_dir_OMNI(date_array[3], data_dir, resolution, file_type)
path = local_dir + filename
df, dict_metadata = read_CDFfile_OMNI(path)
'''
