#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
#matplotlib.rcParams['text.usetex'] = True
#plt.rcParams['font.family'] = 'serif'


'''
###############################################################################
### En este script están todas las funciones para visualizar el campo magnético
### global
################################################################################
'''

def plot_sphere(ax, color='b', alpha=0.7):
    # Crear los parámetros de la esfera
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    # Convertir coordenadas esféricas a coordenadas cartesianas
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))

    # Dibujar la superficie de la esfera
    ax.plot_surface(x, y, z, color=color, alpha=alpha)


def field_line(posit_fl, dens):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x_geo = posit_fl[::dens, 0]
    y_geo = posit_fl[::dens, 1]
    z_geo = posit_fl[::dens, 2]
    ax.scatter(x_geo, y_geo, z_geo)

    plot_sphere(ax)
    ax.set_xlabel('x GEO')
    ax.set_ylabel('y GEO')
    ax.set_zlabel('z GEO')
    ax.set_aspect('equal')

    return ax


def drift_shell(posit, Nposit, dens,dens2, title = ''):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i in range(len(Nposit))[::dens2]:
        field_line = posit[i]

        x_geo = field_line[::dens, 0]
        y_geo = field_line[::dens, 1]
        z_geo = field_line[::dens, 2]
        ax.scatter(x_geo, y_geo, z_geo, color ='dimgrey')

    plot_sphere(ax)
    ax.set_xlabel('x GEO', fontsize=17)
    ax.set_ylabel('y GEO', fontsize=17)
    ax.set_zlabel('z GEO', fontsize=17)
    ax.set_aspect('equal')
    ax.set_title(title)
    ax.tick_params(axis='both', labelsize=16)
    return ax
