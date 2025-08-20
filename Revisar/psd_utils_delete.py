import os, sys
import datetime
import pandas as pd
import time
import numpy as np
import IRBEM

from scipy.interpolate import make_smoothing_spline, UnivariateSpline
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as u
import matplotlib.dates as mdates

from PSD_algorithm_utils import calculate_K_model, calculate_K_integration
from PSD_algorithm_utils import interpolator_alpha_K, interpolate_alpha_K
from PSD_algorithm_utils import calculate_mu, calculate_Emu
from PSD_algorithm_utils import f1_j_alpha, f2_j_alpha, fit_j_alpha
from PSD_algorithm_utils import calculate_j_alphaK, pl_j_energy, exp_f_energy
from PSD_algorithm_utils import flux_to_PSD, calculate_Lstar
from PSD_algorithm_utils import plot_j_alpha_data, plot_flux_psd_energy_data

from create_inputs_utils import get_df_mag_inputs, filter_dates, get_GEO_coordinates
from create_inputs_utils import dict_x_input, dict_magnetic_input_T89

# from create_inputs_utils import get_df_x_inputs



'''
################################################################################
### Definiendo directorios importantes
################################################################################
'''
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

data_omni_dir = home_path + rb_path + data_omni_path


data_rept_dir = home_path + rb_path + data_rept_path





def step0_get_inputs(start_time, dt_hours, df_x_input, fedu, df_mag_input, fedu_metadata, Re):
    print('-------------------------------------------------------------------')
    print('STEP 0 inputs: Get x and magnetic inputs for model.')
    print('-------------------------------------------------------------------')# esto lo entrega para la posición 0 del df
    print('* Initial time analysis: ', start_time)
    print('* Time period of analysis: %s hours' % str(dt_hours))
    print()
    # LAS POSICIONES DE REPT ESTÁN EN GEO, EN KM, TENGO QUE NORMALIZAR POR EL RADIO DE LA TIERRA

    stime = pd.to_datetime(start_time)
    etime = pd.to_datetime(start_time) + pd.Timedelta(hours = dt_hours)

    start_index = df_x_input['Epoch'].searchsorted(stime, side='right')
    end_index = df_x_input['Epoch'].searchsorted(etime, side='right')
#    print(start_index)
#    print(end_index)

    ### Filtramos datos
#    print(len(df_x_input))
    df_x_input_filt = df_x_input.iloc[start_index:end_index]
    list_fedu_filt = fedu[start_index:end_index]

    ###  Recuperamos tiempo y posiciciones (medidas en RE) de los datos REPT
    dates, x_geo, y_geo, z_geo = get_GEO_coordinates(df_x_input_filt, Re)

    ### Iteramos sobre todas las fechas que cumplen las restricciones
    N = len(dates)
    list_x_inputs = []
    list_mag_inputs = []
    for i in range(N):
        ### Creamos diccionario con x_inputs
        x_input = dict_x_input(dates.iloc[i], x_geo[i], y_geo[i], z_geo[i])
        list_x_inputs.append(x_input)

        ### Creamos diccionario con magnetic inputs para T89
        mag_input = dict_magnetic_input_T89(dates.iloc[i], df_mag_input, True)
        list_mag_inputs.append(mag_input)

    #.loc es el índice
    #.iloc es la posición

    return list_x_inputs, list_mag_inputs, list_fedu_filt, fedu_metadata, N


def step5_6_get_f_Emu(js_alphasK, energy_bins, bool, fit_opt, func_opt, alphasK, Esmu, c_unit, flag_plot):

    list_psd = []
    out =[]
    for i in range(len(alphasK)):

        psd = flux_to_PSD(js_to_fit_i, energy_to_fit_i, c_unit)
        psd_unit = psd.unit

        if np.all(np.isnan(js_to_fit_i)):
            list_psd.append(np.nan)
            out.append(np.nan)
            print(psd)
        elif np.sum(bool)< 4:
            print('3 puntos o menos, pasamos')
            list_psd.append(np.nan)
            out.append(np.nan)

        else:
            if fit_opt =='flux':
        #        print('FLUX')
                if flag_plot:
                    ### Ploteamos flujos calculados en alphaK, como función de la energía
                    plot_flux_psd_energy_data(energy_to_fit_i, js_to_fit_i, 'flux', alphasK[i], log=True)

                if func_opt == 'pl':
                    ### Hacemos el ajuste

                    par_opt, _ = curve_fit(pl_j_energy, np.log10(energy_to_fit_i.value), np.log10(js_to_fit_i.value), p0=p_0, bounds = lims)
                    out.append(par_opt)
                    ### Calculamos el valor del flujo para la energía target Emu.
                    log_j_Emu_i = pl_j_energy(np.log10(Esmu[i].value), par_opt[0], par_opt[1])

                elif func_opt == 'spl':
                    #spl = make_smoothing_spline(energy_to_fit.value, js_to_fit.value, lam=0.1)
                    spl = UnivariateSpline(energy_to_fit_i.value, np.log10(js_to_fit_i.value), s=0)
                    out.append(spl)
                    ### Calculamos el valor del flujo para la energía target Emu.
                    log_j_Emu_i = spl(Esmu[i].value)
                else:
                    pass

                ### Transformamos flujo at target Emu a PSD
                j_Emu_i = 10**log_j_Emu_i*js_to_fit_i.unit
                psd_Emu_i = flux_to_PSD(j_Emu_i, Esmu[i], c_unit)
                psd_unit = psd_Emu_i.unit

                if psd_Emu_i.value < 10**(-15):
                    list_psd.append(np.nan)
                else:
                    list_psd.append(psd_Emu_i.value)

            elif fit_opt =='psd':
        #        print('PSD')
                ### Pasamos los flujos calculados para alphaK a PSD
                psd_to_fit_i = flux_to_PSD(js_to_fit_i, energy_to_fit_i, c_unit)
            #    print(psd_to_fit_i)
                if flag_plot:
                    ### Ploteamos flujos calculados en alphaK, como función de la energía
                    plot_flux_psd_energy_data(energy_to_fit_i, psd_to_fit_i, 'psd', alphasK[i], log=True)

                if func_opt == 'exp':

                    ### Hacemos el ajuste
                    p_0 = [10,1]
                    lims = ([0,0], [np.inf, np.inf])
                    par_opt, _ = curve_fit(exp_f_energy, energy_to_fit_i, psd_to_fit_i, p0=p_0, bounds = lims)
                    out.append(par_opt)
                    ### Calculamos el valor de PSD para la energía target Emu
                    psd_Emu_i = exp_f_energy(Esmu[i].value, par_opt[0], par_opt[1])
                    psd_unit = psd_to_fit_i.unit
                if psd_Emu_i < 10**(-15):
                    list_psd.append(np.nan)
                else:
                    list_psd.append(psd_Emu_i)

    list_psd = list_psd
    print('* Target alphasK: ', alphasK)
    print('* Target Emu: ', Esmu)
    print('* PSD at target alphaK and Emu: ', list_psd)
    print()
    return list_psd*psd_unit, out




def checking_step2(Ks_lin_calc, Ks_geom_calc, model_info, calc_opt, Re, spline_lin, spline_geom, Ks, alphasK, flag180=False):
    # Este fragmento de código es para verificar visualmente que la interpolación
    # del step 2 esté funcionando

    fig, ax = plt.subplots()

    [model, x_input, mag_input ]= model_info
    Ks_lin = np.linspace(Ks_lin_calc[-1], Ks_lin_calc[0], 1000)
    alphas_lin = spline_lin(np.log10(Ks_lin.value))

    try:
        Ks_geom = np.linspace(Ks_geom_calc[-1], Ks_geom_calc[0], 1000)
        alphas_geom = spline_geom(np.log10(Ks_geom.value))
    except IndexError:
        Ks_geom = []*Ks_lin.unit
        alphas_geom = []
    ax.plot(np.log10(Ks_lin.value), alphas_lin, color = 'r', label = 'Interp')
    ax.plot(np.log10(Ks_geom.value), alphas_geom, color='r')

    try:
        K_range = np.linspace(Ks_lin_calc[-1], Ks_geom_calc[0], 100)
    except IndexError:
        K_range = np.linspace(Ks_lin_calc[-1], Ks_lin_calc[0], 100)


    K_1 = K_range[5]
    K_2 = K_range[30]
    K_3 = K_range[70]
    alpha_1 = interpolate_alpha_K(K_1, Ks_lin_calc[0], spline_lin, spline_geom, flag180)
    alpha_2 = interpolate_alpha_K(K_2, Ks_lin_calc[0], spline_lin, spline_geom, flag180)
    alpha_3 = interpolate_alpha_K(K_3, Ks_lin_calc[0], spline_lin, spline_geom, flag180)


    if flag180:
        alphas_full = np.linspace(1, 180, 100)*u.degree
        alphas_lin_2 = 180 - alphas_lin
        try:
            alphas_geom_2 = 180 - alphas_geom
        except TypeError:
            alphas_geom_2 = []

        ax.plot(np.log10(Ks_lin.value), alphas_lin_2, color = 'cyan', label = 'Interp 180')
        ax.plot(np.log10(Ks_geom.value), alphas_geom_2, color='cyan')

        ax.plot(np.log10(K_1.value), alpha_1[0].value, 'o', color = 'r')
        ax.plot(np.log10(K_2.value), alpha_2[0].value, 'o', color = 'r')
        ax.plot(np.log10(K_3.value), alpha_3[0].value, 'o', color = 'r')

        ax.plot(np.log10(K_1.value), alpha_1[1].value, 'o', color = 'r')
        ax.plot(np.log10(K_2.value), alpha_2[1].value, 'o', color = 'r')
        ax.plot(np.log10(K_3.value), alpha_3[1].value, 'o', color = 'r')

    else:
        alphas_full = np.linspace(1, 90, 100)*u.degree
        ax.plot(np.log10(K_1.value), alpha_1.value, 'o', color = 'r')
        ax.plot(np.log10(K_2.value), alpha_2.value, 'o', color = 'r')
        ax.plot(np.log10(K_3.value), alpha_3.value, 'o', color = 'r')

    _, Ks_full, _ = calculate_K_model(alphas_full, model, x_input, mag_input, Re, calc_opt)
    ax.plot(np.log10(Ks_full.value), alphas_full.value, '*', color = 'b', label = 'Data')


    for i in range(len(Ks.value)):
        ax.plot(np.log10(Ks[i].value), alphasK[i].value, 'o', color = 'g')

    ax.set_ylabel('alphas')
    ax.set_xlabel('log10(K)')

    ax.set_yscale('log')
    ax.legend()
    plt.show()


def checking_step3(b_mag, alphasK, Esmu, mu):
    E_min = 1*u.MeV # MeV
    E_max = 25*u.MeV # MeV

    Es= np.linspace(E_min, E_max, 1000)
#    print(Es)
    for i in range(len(alphasK.value)):
        mus= calculate_mu(Es, b_mag, alphasK[i])

        fig, ax = plt.subplots()
        ax.plot(Es.value, mus.value, '*', color = 'b')
        ax.set_xlabel('Energy')
        ax.set_ylabel('mu')

        ax.plot(Esmu[i].value, mu, 'o', color = 'r', label = 'Interp')
        plt.show()


def checking_step4(fedu, alpha_bins, energy_bins, fit_flag, parameters, max_values, func_opt, alphasK, js_alpha_K):
    N= len(energy_bins)
    alphas = np.linspace(10,170, 1000)*u.degree
    alphas_rad = np.radians(alphas)

    for i in range(N):
        fig, ax = plt.subplots()
        if fit_flag[i]:
            title_plot = 'Fitted data. Function ' + func_opt
            par_opt = parameters[i]
            max_val = max_values[i]

            if func_opt =='1':
                j = f1_j_alpha(alphas_rad, par_opt[0], par_opt[1], par_opt[2])
            if func_opt =='2':
                j = f2_j_alpha(alphas_rad, par_opt[0], par_opt[1])

            ax.plot(alphas, j*max_val, color = 'r')


        else:
            ax.set_title('No fit')
        ax.set_ylabel('flux [%s]' % str(fedu.unit))
        ax.set_xlabel('Pitch angle [%s]' % str(alpha_bins.unit))

        flux = fedu[:,i]
        ax.plot(alpha_bins, flux, '*', color = 'b', label = energy_bins[i])

        for j in range(len(alphasK)):
            ax.plot(alphasK[j], js_alpha_K[j][i], 'o', color = 'r', label = 'Fitted flux')
        ax.legend()
        ax.set_yscale('log')
        plt.show()
    return


def checking_step5_6(energy_bins, bool,fit_flag, parameters, func_opt, js_alphasK, Esmu, psd_Esmu, c_unit):

    energy = np.linspace(1,20, 100)
    log_energy = np.log10(energy)
    energy_to_plot = energy_bins.value[bool]

    for i in range(len(Esmu)):
        js_to_plot = js_alphasK[i][bool]
        psd_data = flux_to_PSD(js_to_plot, energy_to_plot*energy_bins.unit, c_unit)
        fig, ax = plt.subplots()
        ax.set_yscale('log')
        if np.isnan(parameters[i]):
            print('NAN')
            psd = [np.nan]*len(energy)*psd_Esmu.unit

        else:
            if fit_flag =='flux':

                if func_opt == 'pl':
                    log_j = pl_j_energy(log_energy, parameters[i][0], parameters[i][1])

                elif func_opt == 'spl':
                    energy = np.linspace(energy_to_plot[0], energy_to_plot[-1], 100)
                    log_j = parameters[i](energy)

                else:
                    pass

                ### Transformamos flujo a PSD
                j = 10**log_j*js_to_plot.unit
                psd = flux_to_PSD(j, energy*energy_bins.unit, c_unit)
            elif fit_flag =='psd':
                if func_opt == 'exp':
                    print(parameters[i])
                    psd = exp_f_energy(energy, parameters[i][0], parameters[i][1])

        ax.plot(energy, psd, color = 'r')
        ax.plot(energy_to_plot, psd_data, '*', color = 'b', label = 'Fitted PSD data')
        ax.set_ylabel('PSD [%s]' % str(psd_data.unit))
        ax.set_xlabel('Energy [%s]' % str(energy_bins.unit))

        ax.plot(Esmu[i], psd_Esmu[i], 'o', color = 'r', label = 'Target energy')
        plt.show()


def plot_psd_lstar(lstar_data, psd_data,  start_times, colors, K, mu):
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 7)

    for i in range(len(lstar_data)):
        print(i)
        lstar = np.array(lstar_data[i])
        psd = np.array(psd_data[i])

        mask = psd[1:]/psd[:-1] > 10**(1.5)
        mask2 = np.append(mask, [False])
        psd = psd[~mask2]*psd_data[i].unit
        lstar = lstar[~mask2]

        mask = psd[1:]/psd[:-1] < 10**(-1.5)
        mask2 = np.insert(mask, [0], True)
        psd = psd[~mask2]
        lstar = lstar[~mask2]

        ax.plot(lstar, psd, '-o', markersize= 4, color = colors[i], label = start_times[i])

    title = 'K = ' + str(K) + ',    $\\mu$ = ' + str(mu)
    ax.set_yscale('log')
    ax.set_ylabel('PSD [%s]' % str(psd.unit), fontsize=17)
    ax.set_xlabel('$L^*$', fontsize=17)
    ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
    ax.tick_params(axis='both', labelsize=16)
    ax.set_title(title, fontsize = 20)
    fig.savefig('plot4.pdf')

#    ax.set_xlim([1,10])
    plt.show()
    plt.close(fig)



def plot_psd_time(dates_data, psd_data, start_times, colors, K, mu):

    for i in range(len(psd_data)):
        print(i)
        fig, ax = plt.subplots()
        fig.set_size_inches(8, 7)
        time = dates_data[i]
        psd = psd_data[i]
        title = start_times[i][:10]+ ' ' + str(K) + ' ' + str(mu)
        ax.plot(time, psd, '-o', color = colors[i])

        ax.set_yscale('log')
        ax.set_ylabel('Electron PSD [%s]' % str(psd.unit), fontsize=17)
        ax.set_xlabel('Time', fontsize=17)
#        ax.legend(fancybox=True, shadow=True, ncol=1, fontsize=13)
        ax.tick_params(axis='both', labelsize=16)
        ax.xaxis.set_major_locator(mdates.HourLocator(interval = 1))  # Localizador de días
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))  # Formato de fecha
        ax.set_title(title, fontsize = 20)
        fig.savefig('plot5_'+ str(i)+'.pdf')

#    ax.set_xlim([1,10])
        plt.show()
        plt.close(fig)
