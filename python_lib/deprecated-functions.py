import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from uncertainties import ufloat
import uncertainties
from uncertainties.unumpy import uarray
from scipy.optimize import curve_fit
import os

# Deprecated functions that should be reviewed before being copied to helpers.py


def autoplot(xValues, yValues, xLabel, yLabel, plotLabel="", errorbars=True, plotStyle='ro', errorStyle='g,', yScale='linear', **furtherPlotArgs):
    """Return a subplot object.

    :param errorbars=True: Plots error bars when true.
    :param yScale: e.g. 'log', 'dec'
    """
    xValues = np.array(xValues)
    yValues = np.array(yValues)
    errX = None
    errY = None
    if type(xValues[0]) == uncertainties.Variable or type(xValues[0]) == uncertainties.AffineScalarFunc:
        x = [item.nominal_value for item in xValues]
        errX = [item.std_dev for item in xValues]
    else:
        x = xValues
    if type(yValues[0]) == uncertainties.Variable or type(yValues[0]) == uncertainties.AffineScalarFunc:
        y = [item.nominal_value for item in yValues]
        errY = [item.std_dev for item in yValues]
    else:
        y = yValues

    ax.set_yscale(yScale)
    x_offset = (max(x) - min(x)) * 0.015
    ax.set_xlim(min(x) - x_offset, max(x) + x_offset)
    if yScale != 'log':
        y_offset = (max(y) - min(y)) * 0.015
        ax.set_ylim(min(y) - y_offset, max(y) + y_offset)

    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)

    ax.legend(loc='best')

    if errorbars:
        if errX != None and errY != None:
            plt.errorbar(x, y, xerr=errX, yerr=errY, fmt=errorStyle)
        elif errY != None:
            plt.errorbar(x, y, yerr=errY, fmt=errorStyle)
            print(errY)
        elif errX != None:
            plt.errorbar(x, y, xerr=errX, fmt=errorStyle)
        else:
            raise "Should draw errorbars, but x, y are not ufloats!"
    ax.plot(x, y, plotStyle, label=plotLabel, **furtherPlotArgs)

    fig.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
    return fig



# find peaks
import sys
from numpy import NaN, Inf, arange, isscalar, asarray, array


def peakdet(v, delta, x=None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html

    Returns two arrays

    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.

    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.

    """
    maxtab = []
    mintab = []

    if x is None:
        x = arange(len(v))

    v = asarray(v)

    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')

    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')

    if delta <= 0:
        sys.exit('Input argument delta must be positive')

    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN

    lookformax = True

    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx - delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn + delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return array(maxtab), array(mintab)

# if __name__=="__main__":
#    from matplotlib.pyplot import plot, scatter, show
#    series = [0,0,0,2,0,0,0,-2,0,0,0,2,0,0,0,-2,0]
#    maxtab, mintab = peakdet(series,.3)
#    plot(series)
#    scatter(array(maxtab)[:,0], array(maxtab)[:,1], color='blue')
#    scatter(array(mintab)[:,0], array(mintab)[:,1], color='red')
#    show()


def getPeakVal(peaksmax):
    """gets the values of the peaks for the x and y axes"""
    peakst = []
    for i in range(len(peaksmax)):
        peakst.append(peaksmax[i][0])
    peaksT = []
    for i in range(len(peaksmax)):
        peaksT.append(peaksmax[i][1])
    return peakst, peaksT
