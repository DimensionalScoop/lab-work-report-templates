import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from numpy import *
from scipy.constants import Avogadro, e, electron_mass, hbar, k, mu_0
import scipy.optimize as optimize
import pprint

import data as d
import helpers as hel
import plot_helpers as plot
import table
import simplex



# Amplifier Bandpass Curve
def gaus(x, mu, sigma, scale, noise_level):
    return noise_level + scale * 1 / sqrt(2 * pi * sigma) * exp(-(x - mu)**2 / 2 / sigma**2)

initial_guess = [36e3, 100, 10, 0.1]
params, quality = simplex.optimize(d.amplifier_f, d.amplifier_U, gaus, p0=initial_guess)
fit = lambda x: gaus(x, *params)

U_max = max(d.amplifier_U)
U_half = U_max / sqrt(2)
f_with_max_U = d.amplifier_f[argwhere(d.amplifier_U == U_max)]
half_U_freq_left = optimize.brentq(lambda x: fit(x) - U_half, 34e3, f_with_max_U)
half_U_freq_right = optimize.brentq(lambda x: fit(x) - U_half, f_with_max_U, 38e3)
Q_exp = f_with_max_U / (half_U_freq_right - half_U_freq_left)
print("Bandpass Quality:", Q_exp)

plot.plot(d.amplifier_f, d.amplifier_U,
          "Frequenz in kHz", "Durchgelassene Spannung in V", "amplifier.pdf",
          fit)


# Calculating CHI from measurements
def χ_voltage(sample_voltage, sample_cross_section):
    """Uses bridge's difference in voltage"""

    # 'Exact' Formular doesn't lead to significantly better results
    # f = 35.16e3
    # omega = 2 * pi * f
    # radicant = d.coil_res + omega**2 * (mu_0 * d.coil_n**2 / d.coil_lenght * d.coil_cross_secion)**2
    # return sample_voltage / d.voltage_in * 4 * d.coil_lenght / omega / mu_0 / d.coil_n**2 / sample_cross_section * sqrt(radicant)

    # Approximation Formular
    return sample_voltage / (d.voltage_in) * 4 * d.coil_cross_secion / sample_cross_section


def χ_res(ΔR, sample_cross_section):
    """Uses bridge's resistor matching"""
    return 2 * ΔR / d.res_3 * d.coil_cross_secion / sample_cross_section


χ_results_voltage = []
χ_results_res = []

for element in range(len(d.all_measurements)):

    measurements = d.all_measurements[element]

    voltages = []
    resistances = []
    χ_tmp_voltage = []
    χ_tmp_res = []
    for measurement in measurements:
        voltages.append(measurement.voltages[1])
        resistances.append(measurement.ΔR)

        χ_tmp_voltage.append(χ_voltage(measurement.voltages[1], d.cross_section_areas[element]))
        χ_tmp_res.append(χ_res(measurement.ΔR, d.cross_section_areas[element]))

    χ_results_voltage.append(hel.combine_measurements(χ_tmp_voltage))
    χ_results_res.append(hel.combine_measurements(χ_tmp_res))

    table.petit_table([χ_tmp_voltage, array(voltages) * 1e6, χ_tmp_res, resistances], "data_element_" + str(element), [4, 4, 4, 6])

χ_mean = hel.cutErrors((array(χ_results_voltage) + array(χ_results_res)) * 0.5)


# Calculating CHI with HUND
def calc_g_J(S, L, J):
    return (3 * J * (J + 1) + (S * (S + 1) - L * (L + 1))) / (2 * J * (J + 1))


def χ_hund(S, L, J, N):
    T = 288.15  # kelvin at room temperature (15°C)
    mu_B = 0.5 * e / electron_mass * hbar
    g_J = calc_g_J(S, L, J)
    return mu_0 / k * mu_B * N * mu_B * g_J**2 * J * (J + 1) / (3 * T)

element_names = ["$\\text{Dy}^{3+}$", "$\\text{Nd}^{3+}$", "$\\text{Gd}^{3+}$"]
hund_params = [[2.5, 5, 7.5], [1.5, 6, 4.5], [3.5, 0, 3.5]]  # [Spin, Magnetzahl, Gesamtdrehimpuls, Teilchendichte]

χ_reference = [χ_hund(*hund_params[i], d.amounts[i]) for i in range(len(hund_params))]
χ_U_deviation = [hel.abweichung(mean, ref) for mean, ref in zip(χ_results_voltage, χ_reference)]
χ_U_deviation = hel.cutErrors(χ_U_deviation)
χ_R_deviation = [hel.abweichung(mean, ref) for mean, ref in zip(χ_results_res, χ_reference)]
χ_R_deviation = hel.cutErrors(χ_R_deviation)
χ_deviation = [hel.abweichung(mean, ref) for mean, ref in zip(χ_mean, χ_reference)]

# Debugging
print("Get those values to 0:")
print(sum(χ_deviation), sum(χ_U_deviation), sum(χ_R_deviation))

# Generate tables
g_J = [calc_g_J(*par) for par in hund_params]
electron_conf = [[*param, g] for param, g in zip(hund_params, g_J)]

hund_table = [χ_results_voltage, χ_results_res, χ_mean, χ_reference, χ_U_deviation, χ_R_deviation, χ_deviation]
table.write_table(hund_table, "../tables/chi-results.tex", figures=[1, 1, 5, 5, 3, 3, 3],
                  row_names=element_names)

table.write_table(electron_conf, "../tables/hund_params.tex", figures=[2, 2, 2, 4],
                  row_names=['$\\text{Spinsumme}$', '$\\text{Magnetzahlensumme}$', '$\\text{Gesamtdrehimpuls}$', '$\\text{Lande-Faktor}$'])  # , '$\\text{Teilchendichte in} 1/m^3$'])

# Parameter der Proben
table.petit_table(transpose([d.amounts, d.lenghts_total, d.lenghts_inserted, d.cross_section_areas]), "proben_parameter", [3, 3, 3], dont_align=True,
                  row_names=['$\\text{Seltene-Erden-Atome pro m}^3$', '$\\text{Probenlänge in m}$', '$\\text{Eingeführte Probenlänge in m}$', '$\\text{Effektiver Querschnitt in m}^2$'])



# Plot results
plt.clf()
x = [1, 2, 3]
y1, err1 = plot.extract_error(χ_results_voltage)
y2, err2 = plot.extract_error(χ_results_res)
y3 = χ_reference
y4 = χ_mean
err4 = err1 + err2
labels = ["$\mathrm{Dy}_2\mathrm{O}_3$", "$\mathrm{Gd}_2\mathrm{O}_3$", "$\mathrm{Nd}_2\mathrm{O}_3$"]

plt.errorbar(x, y1, yerr=err1, fmt='b.', label="Spannungsmessung")
plt.errorbar(x, y2, yerr=err2, fmt='r.', label="Wiederstandsmessung")
plt.plot(x, y3, 'gx', label="Hund'sche Regel")
plt.plot(x, y4, 'c_', label="Mittelwert beider Messungen")


plt.xlim(0, 4)
plt.ylim(*plot.autolimits(append(append(y1, y2), y3), err=append(err1, err2)))

plt.xticks(x, labels)

plt.xlabel("Probe")
plt.ylabel("$\chi$")
plt.legend(loc='best')

plt.savefig(plot.plot_path() + "chi.pdf")
