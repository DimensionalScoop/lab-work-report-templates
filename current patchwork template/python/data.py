import copy

import numpy as np
import uncertainties
from numpy import array
from scipy.constants import Avogadro, u
from uncertainties import unumpy as unp
from uncertainties import ufloat

import functions as f
from table import petit_table

# Band Pass and Amplifier
amplifier_f = np.array([
    4026, 3940, 3877, 3830, 3810, 3747, 3719, 3671, 3596, 3615, 3567,
    3508, 3548, 3519, 3490, 3500, 3453, 3390, 3258, 3125, 3503, 3506,
    3512, 3514.5]) * 1e1

amplifier_U = np.array([
    39.5, 47, 55, 63, 68, 81.5, 99, 130, 265, 210, 430, 750, 690, 990, 450, 570, 240,
    130, 65, 43, 600, 660, 770, 845]) * 1e-3


# Cross Sections
density_Dy = 7.8e3
density_Gd = 7.4e3
density_Nd = 7.24e3
densities = array([density_Dy, density_Nd, density_Gd])

masses_total = array([185e-4, 9e-3, 1408e-5])
lenghts_total = array([175e-3, 172e-3, 176e-3])
lenghts_inserted = lenghts_total - array([2e-2, 25e-3, 196e-4])
lenght_rel_difference = lenghts_inserted / lenghts_total

cross_section_areas = masses_total / (lenghts_total * densities) * lenght_rel_difference


#masses_adjusted = densities * (lenghts * cross_section_areas)
molar_masses = array([373, 362.5, 336.48]) * 1e-3
amounts = 2 * densities / molar_masses * Avogadro  # number of particles per volume in a sample

# Bridge Parameters
coil_cross_secion = 86.6e-6
coil_n = 250
coil_lenght = 135e-3
coil_res = 0.7  # ohm
voltage_in = 1
res_3 = 998   # ohm



# Bridge Measurements


class measurement():

    def __init__(self, voltages, res_before, res_after):
        """Converts to the correct SI units automatically"""
        self.voltages = array(voltages) * 1e-3 * 1e-2  # amplifier shifts for 1e2
        self.res_before = res_before * 5e-3
        self.res_after = res_after * 5e-3
        self.ΔR = self.res_before - self.res_after

dy_mesurements = [
    measurement([2.45, 59.5, 3.75], 587.6, 191.2),
    measurement([2.6, 55, 3.95], 557.6, 178.4),
    measurement([2.5, 53.5, 4.2], 569.6, 199.6)
]

nd_mesurements = [
    measurement([2.3, 4.15, 2.5], 551.6, 538.2),
    # discarded one measurement, as written in data log
    measurement([2.15, 4.6, 2.2], 549.2, 525.7),
    measurement([2.18, 3.4, 2.2], 543.6, 525.6)
]

gd_mesurements = [
    measurement([2.2, 24.5, 2.9], 568.5, 395.5),
    measurement([2.2, 23, 2.95], 560, 400),
    measurement([2.2, 24, 2.9], 564.2, 402.5)
]

all_measurements = [dy_mesurements, nd_mesurements, gd_mesurements]


# Reference Values
χ_reference = [0, 0, 0]
