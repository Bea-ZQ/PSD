###############################################################################
### En este script estoy probando el modelo Tsyganeno 89 y todas sus rutinas
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

'''
################################################################################
# Probando make_lstar float
################################################################################
date= '2013-05-30T00:00:10.254'
pd_date = pd.to_datetime(date)
x1 = 5.16788058
x2 = -0.96177772
x3 = 0.89858463

x_input = {'x1':x1, 'x2':x2, 'x3':x3, 'dateTime':date}
mag_input = {'Kp':40}

out_dict = model_t89.make_lstar(x_input, mag_input)

out_dict_gfm = model_t89.get_field_multi(x_input, mag_input)
b_mag = out_dict_gfm['Bl']

out_dict_fmeq = model_t89.find_magequator(x_input, mag_input)
b_min = out_dict_fmeq['bmin']

mlt = model_t96.get_mlt(x_input)

'''


'''
################################################################################
# Probando get MLT
################################################################################
#                       Epoch       MLAT        MLT      B_calc         I         L    L_star   Position_x1   Position_x2  Position_x3
#0    2013-05-30 00:00:10.254  13.755815  23.419799  195.272041  1.764406  6.124262  5.521910  32924.567151  -6127.485885  5724.882676
#1    2013-05-30 00:00:21.362  13.757659  23.422461  195.477909  1.763058  6.121805  5.520014  32911.548768  -6129.907717  5723.246079
#7775 2013-05-30 23:59:33.360        NaN  21.202239  204.375576  2.356991  6.293915  5.676600  24582.647235 -23198.851964  4550.357392
#7776 2013-05-30 23:59:44.467        NaN  21.204867  204.198975  2.360597  6.296957  5.678981  24587.548319 -23210.696741  4554.310374

#array([[ 5.16788058, -0.96177772,  0.89858463],
#       [ 5.16583719, -0.96215786,  0.89832775],
#       [ 3.85852256, -3.64132035,  0.7142297 ],
#       [ 3.85929184, -3.64317952,  0.71485016]])


date1= '2013-05-30T00:00:10.254'
x1_input = {'x1':5.16788058, 'x2':-0.96177772, 'x3':0.89858463, 'dateTime':date1}
mlt1 = model_t96.get_mlt(x1_input)

date2= '2013-05-30T00:00:21.362'
x2_input = {'x1':5.16583719, 'x2':-0.96215786, 'x3':0.89832775, 'dateTime':date2}
mlt2 = model_t96.get_mlt(x2_input)
'''



'''
################################################################################
# Probando trace_field_line: LISTO
################################################################################

date= '2013-05-30T00:00:10.254'
pd_date = pd.to_datetime(date)
x1 = 5.16788058
x2 = -0.96177772
x3 = 0.89858463
x_input = {'x1':x1, 'x2':x2, 'x3':x3, 'dateTime':date}
mag_input = {'Kp':40}
R0 = 3
dens = 2
out_dict = model_t89.trace_field_line(x_input, mag_input, R0)
posit = out_dict['POSIT']
Nposit = out_dict['Nposit']

# Now plot the field lines
ax = plot_field_line(posit, dens)
ax.scatter(x1, x2, x3, color='red', s=100, label = 'input pos')
ax.legend()
plt.show()

blocal = out_dict['blocal']
bmin = out_dict['bmin']

out_dict_fmeq = model_t89.find_magequator(x_input, mag_input)
b_min_fmeq = out_dict_fmeq['bmin']

bs_mag = []
for i in range(len(posit)):
    x_geo_i = posit[i][0]
    y_geo_i = posit[i][1]
    z_geo_i = posit[i][2]

    x_input_i = {'x1':x_geo_i, 'x2':y_geo_i, 'x3':z_geo_i, 'dateTime':date}
    out_dict_i = model_t89.get_field_multi(x_input_i, mag_input)
    b_mag = out_dict_i['Bl'][0]
    bs_mag.append(b_mag)
'''

'''
################################################################################
# Probando find_mirror_points: LISTO
################################################################################

date= '2013-05-30T00:00:10.254'
pd_date = pd.to_datetime(date)
x1 = 5.16788058
x2 = -0.96177772
x3 = 0.89858463
x_input = {'x1':x1, 'x2':x2, 'x3':x3, 'dateTime':date}
mag_input = {'Kp':40}
alpha = 50
R0 = 1
dens = 5
out_dict_fl = model_t89.trace_field_line(x_input, mag_input, R0)
posit_fl = out_dict_fl['POSIT']

out_dict = model_t89.find_mirror_point(x_input, mag_input, alpha)
posit = out_dict['POSIT']
blocal = out_dict['blocal']
bmin = out_dict['bmin']

x_geo = posit[0]
y_geo = posit[1]
z_geo = posit[2]

out_dict_bl = model_t89.get_field_multi(x_input, mag_input)
x_input_bmirr = {'x1':x_geo, 'x2':y_geo, 'x3':z_geo, 'dateTime':date}
out_dict_bmirr = model_t89.get_field_multi(x_input_bmirr, mag_input)


ax = plot_field_line(posit_fl, dens)
ax.scatter(x1, x2, x3, color='red', s=50, label = 'input pos')
ax.scatter(x_geo, y_geo, z_geo, color='black', s=50, label = 'mirror pos')
ax.legend()
plt.show()
'''
'''
################################################################################
# Probando mirror_point_altitude
################################################################################

date= '2013-05-30T00:00:10.254'
pd_date = pd.to_datetime(date)
x1 = 5.16788058
x2 = -0.96177772
x3 = 0.89858463
x_input = {'x1':x1, 'x2':x2, 'x3':x3, 'dateTime':date}
mag_input = {'Kp':40}

R0 = 1
alpha = 90
dens = 5
out_dict_fl = model_t89.trace_field_line(x_input, mag_input, R0)
posit_fl = out_dict_fl['POSIT']

out_dict_mi = model_t89.find_mirror_point(x_input, mag_input, alpha)
posit_mi = out_dict_mi['POSIT']
b_mi = out_dict_mi['blocal']
x_geo_mi = posit_mi[0]
y_geo_mi = posit_mi[1]
z_geo_mi = posit_mi[2]

altura, posit_mi2 = model_t89.mirror_point_altitude(x_input, mag_input, R0)

opp_hem = 2
out_dict_foot = model_t89.find_foot_point(x_input, mag_input, altura, opp_hem)

foot_gdz = out_dict_foot['XFOOT']
bmag_foot = out_dict_foot['BFOOTMAG']

foot_geo = IRBEM.Coords().transform(time=pd_date, pos = foot_gdz, sysaxesIn='GDZ', sysaxesOut='GEO')




ax = plot_field_line(posit_fl, dens)
ax.scatter(x1, x2, x3, color='red', s=50, label = 'input pos')
ax.scatter(x_geo_mi, y_geo_mi, z_geo_mi, color='black', s=50, label = 'mirror pos')
ax.scatter(foot_geo[0][0], foot_geo[0][1], foot_geo[0][2], color='green', s=50, label= 'foot pos')
ax.legend()
plt.show()

'''
'''
################################################################################
# Probando get_field_multi: LISTO
################################################################################

date_0= '2013-05-30T00:00:10.254'
pd_date_0= pd.to_datetime(date_0)

date_1= '2013-05-30T23:59:44.467'
pd_date_1 = pd.to_datetime(date_1)

dates = np.array([pd_date_0, pd_date_1])

x1_0 = 5.16788058
x2_0 = -0.96177772
x3_0 = 0.89858463

x1_1 = 3.85929184
x2_1 = -3.64317952
x3_1 = 0.71485016

x1 = [x1_0, x1_1]
x2 = [x2_0, x2_1]
x3 = [x3_0, x3_1]

x_input = {'x1':x1, 'x2':x2, 'x3':x3, 'dateTime':dates}

kps = [40,50]
mag_input = {'Kp':kps}

out_dict = model_t89.get_field_multi(x_input, mag_input)

bx_geo = out_dict['BxGEO']
by_geo = out_dict['ByGEO']
bz_geo = out_dict['BzGEO']
b_mag = out_dict['Bl']
'''

'''
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
mag_input = {'Kp':kps}

out_dict = model_t89.find_magequator(x_input, mag_input)

b_min = out_dict['bmin']
x_geo = out_dict['XGEO'][0]
y_geo = out_dict['XGEO'][1]
z_geo = out_dict['XGEO'][2]

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


'''
################################################################################
# Probando drift_shell: Listo
################################################################################


date= '2013-05-30T00:00:10.254'
pd_date= pd.to_datetime(date)

x1 = 5.16788058
x2 = -0.96177772
x3 = 0.89858463
x_input = {'x1':x1, 'x2':x2, 'x3':x3, 'dateTime':date}

kps = 40
mag_input = {'Kp':kps}

dens = 5

out_dict = model_t89.drift_shell(x_input, mag_input)
posit = out_dict['POSIT']
Nposit = out_dict['Nposit']
blocal = out_dict['blocal']
bmin = out_dict['bmin']


ax = plot_drift_shell(posit, Nposit, dens, 6)
ax.scatter(x1, x2, x3, color='red', s=50, label = 'input pos')
ax.legend()
plt.show()


out_dict_fmeq = model_t89.find_magequator(x_input, mag_input)
b_min_fmeq = out_dict_fmeq['bmin']

bs_mag = []
field_line = posit[8]
for i in range(len(field_line)):
    x_geo_i = field_line[i][0]
    y_geo_i = field_line[i][1]
    z_geo_i = field_line[i][2]

    x_input_i = {'x1':x_geo_i, 'x2':y_geo_i, 'x3':z_geo_i, 'dateTime':date}
    out_dict_i = model_t89.get_field_multi(x_input_i, mag_input)
    b_mag = out_dict_i['Bl'][0]
    bs_mag.append(b_mag)

'''



'''
################################################################################
# Probando bounce_period: LISTO
################################################################################

date= '2013-05-30T00:00:10.254'
pd_date= pd.to_datetime(date)

x1 = 5.16788058
x2 = -0.96177772
x3 = 0.89858463
x_input = {'x1':x1, 'x2':x2, 'x3':x3, 'dateTime':date}

kps = 40
mag_input = {'Kp':kps}

alpha = 90
R0=1
E = np.arange(200, 1000, 100)

bounce_periods = model_t89.bounce_period(x_input, mag_input, E, Erest=511, R0=1, alpha=alpha, interpNum=100000)

###

model = IRBEM.MagFields(options = [0,0,0,0,0], verbose = True)

date_2= '2015-02-02T22:00:00'
pd_date= pd.to_datetime(date_2)

x1 = 651
x2 = 65
x3 = 15.9
x_input_2 = {'x1':x1, 'x2':x2, 'x3':x3, 'dateTime':date_2}

kps = 40
mag_input_2 = {'Kp':kps}

E = np.arange(200, 1000)

bounce_periods_2 = model.bounce_period(x_input_2, mag_input_2, E)
'''


################################################################################
# Probando drift_bounce_orbit: LISTO
################################################################################

date= '2013-05-30T00:00:10.254'
pd_date= pd.to_datetime(date)

x1 = 5.16788058
x2 = -0.96177772
x3 = 0.89858463
x_input = {'x1':x1, 'x2':x2, 'x3':x3, 'dateTime':date}

kps = 40
mag_input = {'Kp':kps}

alpha = 60
R0 = 1
out_dict = model_t89.drift_bounce_orbit(x_input, mag_input, alpha, R0)
posit = out_dict['POSIT']
Nposit = out_dict['Nposit']
blocal = out_dict['blocal']
bmin = out_dict['bmin']
bmirr = out_dict['bmirr']

b_min_me = []
flag= True
for i in range(len(Nposit)):
    x_input_me = {'x1':posit[i][0][0], 'x2':posit[i][0][1], 'x3':posit[i][0][2], 'dateTime':date}
    out_dict_me = model_t89.find_magequator(x_input_me, mag_input)
    b_min_me.append(out_dict_me['bmin'])
    if flag:
        x_geo_me = out_dict_me['XGEO'][0]
        y_geo_me = out_dict_me['XGEO'][1]
        z_geo_me = out_dict_me['XGEO'][2]
        flag = False

b_mirr_mi = []
flag= True
for i in range(len(Nposit)):
    x_input_mi = {'x1':posit[i][0][0], 'x2':posit[i][0][1], 'x3':posit[i][0][2], 'dateTime':date}
    out_dict_mi = model_t89.find_mirror_point(x_input_mi, mag_input, alpha)
    b_mirr_mi.append(out_dict_mi['bmin'])
    if flag:
        posit_mi = out_dict_mi['POSIT']
        x_geo_mi = posit_mi[0]
        y_geo_mi = posit_mi[1]
        z_geo_mi = posit_mi[2]
        flag = False


bs_mag = []
field_line = posit[0]
for i in range(len(field_line)):
    x_geo_i = field_line[i][0]
    y_geo_i = field_line[i][1]
    z_geo_i = field_line[i][2]

    x_input_i = {'x1':x_geo_i, 'x2':y_geo_i, 'x3':z_geo_i, 'dateTime':date}
    out_dict_i = model_t89.get_field_multi(x_input_i, mag_input)
    b_mag = out_dict_i['Bl'][0]
    bs_mag.append(b_mag)


dens = 5
ax = plot_drift_shell(posit, Nposit, dens, 3)
ax.scatter(x1, x2, x3, color='red', s=50, label = 'input pos')
ax.scatter(x_geo_me, y_geo_me, z_geo_me, color='blue', s=50, label = 'm_eq pos')
ax.scatter(x_geo_mi, y_geo_mi, z_geo_mi, color='green', s=50, label = 'mirror pos')
ax.legend()
plt.show()
