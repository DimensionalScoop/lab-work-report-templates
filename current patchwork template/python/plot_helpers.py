import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import uncertainties
import uncertainties.unumpy as unp
from scipy.constants import C2K, K2C
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import uarray


def extract_error(data):
    if(isinstance(data[0], uncertainties.UFloat)):
        error = unp.std_devs(data)
        nominal = unp.nominal_values(data)
    else:
        nominal = data
        error = None
    return nominal, error


def autolimits(data, err=None):
    min_lim = min(data)
    max_lim = max(data)
    offset = (max(data) - min(data)) * 0.025
    if err is not None:
        offset += max(err) / 2
    return [min_lim - offset, max_lim + offset]


def plot_path():
    """Returns the correct path to the plots folder in M&Ts praktikum"""
    path, folder = os.path.split(os.getcwd())
    if folder[0] == "V":
        return 'plots/'
    else:
        return '../plots/'


def plot(x_messung, y_messung, xlabel, ylabel, filename, theorie):
    """Plottet diskrete Messwerte gegen eine kontinuierliche Messkurve

    Args:
        x_messung (uarray)
        y_messung (uarray)
        theorie (func(x)): Theoriefunktion, die x-Werte annimmt und y-Werte ausspuckt
        xlabel (string)
        ylabel (string)
        filename (string)

    Returns:
        TYPE: None
    """
    plt.clf()

    x_messung, x_error = extract_error(x_messung)
    y_messung, y_error = extract_error(y_messung)

    x_limit = autolimits(x_messung, err=x_error)
    x_flow = np.linspace(*x_limit, num=1000)
    y_messung = y_messung


    if theorie is not None:
        plt.plot(x_flow, theorie(x_flow), 'g-', label="Fit")
    if x_error is not None:
        plt.errorbar(x_messung, y_messung, xerr=x_error, yerr=y_error, fmt='r,', label="Fehler")

    plt.plot(x_messung, y_messung, 'r.', label="Messwerte")

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.legend(loc='best')


    plt.grid()
    plt.savefig(plot_path() + filename)
