###############################################################################
### En este script estoy probando el modelo Tsyganeno 96 y todas sus rutinas
################################################################################

import os, sys
import datetime
import pandas as pd
import numpy as np
import IRBEM
from plot_utils import plot_field_line, plot_drift_shell
import matplotlib.pyplot as plt

################################################################################
###                        Objetos clase MagFields                           ###
################################################################################
options = [0,0,0,0,0]
options = [1,0,0,0,0]
options = [2,0,0,0,0]
kext_t89 = 'T89'
kext_t96 = 'T96'

model_t89 = IRBEM.MagFields(options=options, kext=kext_t89, verbose = True, sysaxes=1)
model_t96 = IRBEM.MagFields(options=options, kext=kext_t96, verbose = True,  sysaxes=1)

################################################################################
# Probando find_magequator: LISTO
################################################################################


date= '2013-05-30T00:00:10.254'
pd_date= pd.to_datetime(date)

x1 = 5.16788058
x2 = -0.96177772
x3 = 0.89858463

x_input = {'x1':x1, 'x2':x2, 'x3':x3, 'dateTime':date}

kps = 40
dst= -10.0
pdyn = 1.3
by = 5.0
bz = -3.0
mag_input = {'Dst': dst, 'Pdyn': pdyn, 'ByIMF':by, 'BzIMF':bz}

out_dict = model_t96.find_magequator(x_input, mag_input)

#b_min = out_dict['bmin']
#x_geo = out_dict['XGEO'][0]
#y_geo = out_dict['XGEO'][1]
#z_geo = out_dict['XGEO'][2]

'''
R0 = 1
dens = 5
out_dict_fl = model_t89.trace_field_line(x_input, mag_input, R0)
posit_fl = out_dict_fl['POSIT']


ax = plot_field_line(posit_fl, dens)
ax.scatter(x1, x2, x3, color='red', s=50, label = 'input pos')
ax.scatter(x_geo, y_geo, z_geo, color='black', s=50, label = 'mag equator pos')
ax.legend()
plt.show()

bs=[]
for i in range(len(posit_fl)):
    x_geo_i = posit_fl[i][0]
    y_geo_i = posit_fl[i][1]
    z_geo_i = posit_fl[i][2]

    x_input_i = {'x1':x_geo_i, 'x2':y_geo_i, 'x3':z_geo_i, 'dateTime':date}
    out_dict = model_t89.get_field_multi(x_input_i, mag_input)
    b_mag = out_dict['Bl']
    bs.append(b_mag)
'''


'''
################################################################################
# Probando find_foot_point: LISTO
################################################################################


date= '2013-05-30T00:00:10.254'
pd_date= pd.to_datetime(date)

x1 = 5.16788058
x2 = -0.96177772
x3 = 0.89858463

x_input = {'x1':x1, 'x2':x2, 'x3':x3, 'dateTime':date}

kps = 40
mag_input = {'Kp':kps}

same_hem =0
oppo_hem = 2
alt =1000

out_dict_1 = model_t89.find_foot_point(x_input, mag_input, alt, same_hem)
out_dict_2 = model_t89.find_foot_point(x_input, mag_input, alt, oppo_hem)

foot_gdz_1 = out_dict_1['XFOOT']
B_foot_1 = out_dict_1['BFOOT']
Bmag_foot_1 = out_dict_1['BFOOTMAG']

foot_gdz_2 = out_dict_2['XFOOT']
B_foot_2 = out_dict_2['BFOOT']
Bmag_foot_2 = out_dict_2['BFOOTMAG']

foot_geo_1 = IRBEM.Coords().transform(time=pd_date, pos = foot_gdz_1, sysaxesIn='GDZ', sysaxesOut='GEO')
foot_geo_2 = IRBEM.Coords().transform(time=pd_date, pos = foot_gdz_2, sysaxesIn='GDZ', sysaxesOut='GEO')



R0 = 1
dens = 5
out_dict_fl = model_t89.trace_field_line(x_input, mag_input, R0)
posit_fl = out_dict_fl['POSIT']

ax = plot_field_line(posit_fl, dens)
ax.scatter(x1, x2, x3, color='red', s=50, label = 'input pos')

ax.scatter(foot_geo_1[0][0], foot_geo_1[0][1], foot_geo_1[0][2], color='black', s=50, label = 'foot points')
ax.scatter(foot_geo_2[0][0], foot_geo_2[0][1], foot_geo_2[0][2], color='black', s=50)

ax.legend()
plt.show()

x_input_b1 = {'x1':foot_geo_1[0][0], 'x2':foot_geo_1[0][1], 'x3':foot_geo_1[0][2], 'dateTime':date}
x_input_b2 = {'x1':foot_geo_2[0][0], 'x2':foot_geo_2[0][1], 'x3':foot_geo_2[0][2], 'dateTime':date}

out_dict_b1 = model_t89.get_field_multi(x_input_b1, mag_input)
out_dict_b2 = model_t89.get_field_multi(x_input_b2, mag_input)
'''
