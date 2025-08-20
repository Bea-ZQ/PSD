###############################################################################
### En este script están todas las funciones que me permiten crear los inputs
### para los modelos de campo magnético
################################################################################

import os, sys
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def plot_omni(x, y, labely):
    fig, ax = plt.subplots()
    ax.plot(x, y)

    # Formatear el eje x para mostrar solo los días
    ax.xaxis.set_major_locator(mdates.DayLocator(interval = 5))  # Localizador de días
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))  # Formato de fecha

    # Rotar las etiquetas del eje x para que se vean mejor
    plt.xticks(rotation=45)

    # Agregar etiquetas y título
    ax.set_xlabel('Date')
    ax.set_ylabel(labely)
    ax.set_title('Magnetic Inputs Over Time')

    # Mostrar la gráfica
    plt.tight_layout()
    plt.show()


def plot_orbit(coordinates_A, coordinates_B, label, lim):

    def plot(x_A, y_A, x_B, y_B, xlabel, ylabel, xlim, ylim):
        fig, ax = plt.subplots()
        ax.plot(x_A, y_A, label = 'RBSP A', color = 'red')
        ax.plot(x_B, y_B, label = 'RBSP B', color = 'blue')
        tierra = plt.Circle((0, 0), 1, color='gray', fill=True)

        # Agregar etiquetas y título
        ax.set_xlim([-1*xlim, xlim])  # Ajustar los límites de los ejes
        ax.set_ylim([-1*ylim, ylim])

        ax.set_xlabel(xlabel, fontsize=15)
        ax.set_ylabel(ylabel, fontsize=15)
        ax.set_title('RBSP ORBIT')
        ax.add_artist(tierra)
        ax.set_aspect('equal')  # Mantener la proporción 1:1
        ax.legend()
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)
        plt.tight_layout()
        plt.grid(True)
        plt.show()
    x_A = coordinates_A[:,0]
    y_A = coordinates_A[:,1]
    z_A = coordinates_A[:,2]

    x_B = coordinates_B[:,0]
    y_B = coordinates_B[:,1]
    z_B = coordinates_B[:,2]

    plot(x_A, y_A, x_B, y_B, 'X '+ label, 'Y '+ label, lim, lim)
    plot(y_A, z_A, y_B, z_B, 'Y '+ label, 'Z '+ label, lim, lim/2)
    plot(x_A, z_A, x_B, z_B, 'X '+ label, 'Z '+ label, lim, lim/2)



def plot_orbit_3D(x, y, z, sys):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(x, y, z, color = 'b')
    ax.scatter(x[0], y[0], z[0], color = 'r', label='init pos')
    print('X: ', x[0])
    print('Y: ', y[0])
    print('Z: ', z[0])

    plot_sphere(ax)
    ax.set_xlabel('x '+sys)
    ax.set_ylabel('y '+sys)
    ax.set_zlabel('z '+sys)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.set_aspect('equal')
    ax.legend()
    return
