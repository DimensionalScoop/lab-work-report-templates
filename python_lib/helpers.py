import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import uncertainties
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import uarray

from python_lib.simplex import *  # allow simplex fitting by importing only helpers

maxfev = 1000000


def autofit(x, y, fitFunction, p0=None):
    """Returns params of the curvefit as ufloat. y can be an array of ufloats."""
    if isinstance(y[0], uncertainties.UFloat):
        ny = [i.nominal_value for i in y]
        dy = [i.std_dev for i in y]
        params, covariance = curve_fit(fitFunction, x, ny, sigma=dy, absolute_sigma=True,
                                       p0=p0, maxfev=maxfev)
    else:
        params, covariance = curve_fit(fitFunction, x, y, p0=p0, maxfev=maxfev)
    errors = np.sqrt(np.diag(covariance))
    return uarray(params, errors)


def combine_measurements(values):
    """Combines a np.array of measurements into one ufloat by taking the mean."""
    return ufloat(mean(values), stdDevOfMean(values))


def mean(values):
    """Return the mean of values"""
    values = np.array(values)
    return sum(values) / len(values)


def stdDev(values):
    """Return estimated standard deviation"""
    values = np.array(values)
    b = 0
    m = mean(values)
    for x in values:
        b += (x - m) ** 2
    return np.sqrt(1 / (len(values) - 1) * b)


def stdDevOfMean(values):
    """Return estimated standard deviation of the mean (the important one!)"""
    return stdDev(values) / np.sqrt(len(values))


def estimate_sigmas(values, ableseunsicherheit):
    """Generates std deviations for analogue instruments. Returns a ufloatarray."""
    nominal = values
    magnitude = np.floor(np.log10(nominal))
    error = [ableseunsicherheit * 10**mag for mag in magnitude]

    return uarray(nominal, error)


def estimate_sigmas_only(values, ableseunsicherheit):
    """Generates std deviations for analogue instruments. Returns only an array with the errors."""
    nominal = values
    magnitude = np.floor(np.log10(nominal))
    error = [ableseunsicherheit * 10**mag for mag in magnitude]

    return error

def dev_from_theo(value, lit):
    """Returns deviation of an experimental value from a literature value."""
    return (lit - value) / lit * 100
