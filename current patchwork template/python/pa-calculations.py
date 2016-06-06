import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix, ufloat

#-------------------------Ub-180V---------------------------------------------
Ua11 = ufloat(0.000023, 0.000001)
Ua12 = ufloat(0.000024, 0.000001)
Ua13 = ufloat(0.000025, 0.000001)
Ua21 = ufloat(0.000025, 0.000001)
Ua22 = ufloat(0.000027, 0.000001)
Ua23 = ufloat(0.000027, 0.000001)
Ua31 = ufloat(0.000027, 0.000001)
Ua32 = ufloat(0.000026, 0.000001)
Ua33 = ufloat(0.000028, 0.000001)
Ua34 = ufloat(0.000026, 0.000001)

UBr11 = ufloat(0.000340, 0.000001)
UBr12 = ufloat(0.000360, 0.000001)
UBr13 = ufloat(0.000340, 0.000001)
UBr21 = ufloat(0.000730, 0.000001)
UBr22 = ufloat(0.000720, 0.000001)
UBr23 = ufloat(0.000750, 0.000001)
UBr31 = ufloat(0.000060, 0.000001)
UBr32 = ufloat(0.000080, 0.000001)
UBr33 = ufloat(0.000088, 0.000001)
UBr34 = ufloat(0.000083, 0.000001)

dR11 = ufloat((0.534 - 0.352) * 5, 0.001)
dR12 = ufloat((0.522 - 0.333) * 5, 0.001)
dR13 = ufloat((0.525 - 0.350) * 5, 0.001)
dR21 = ufloat((0.529 - 0.142) * 5, 0.001)
dR22 = ufloat((0.527 - 0.136) * 5, 0.001)
dR23 = ufloat((0.534 - 0.140) * 5, 0.001)
dR31 = ufloat((0.523 - 0.500) * 5, 0.001)
dR32 = ufloat((0.532 - 0.494) * 5, 0.001)
dR33 = ufloat((0.525 - 0.487) * 5, 0.001)
dR34 = ufloat((0.529 - 0.494) * 5, 0.001)
# print('rd')
# print(0.352*5)
# print(0.333*5)
# print(0.350*5)
# print(0.142*5)
# print(0.136*5)
# print(0.140*5)
# print(0.500*5)
# print(0.494*5)
# print(0.487*5)
# print(0.494*5)

s1 = 3.5
s2 = 2.5
s3 = 1.5
j1 = 3.5
j2 = 7.5
j3 = 4.5
l1 = 0
l2 = 5
l3 = 6
M1 = 362 * 1.660539040 * 10**(-27)
M2 = 373 * 1.660539040 * 10**(-27)
M3 = 336 * 1.660539040 * 10**(-27)



R3 = ufloat(998, 1)

vm = ufloat(35000, 10)
vp = ufloat(35400, 10)
v0 = ufloat(35200, 10)

omega = 2 * np.pi * v0

Na = 6.022140857 * 10**23  # mol-1
muh0 = 4 * np.pi * 10**(-7)
hq = 6.626070040 * 10**(-34) / (2 * np.pi)
e0m0 = -1.758820024 * 10**11
k = 1.38064852 * 10**(-23)
T = 298.15
muhB = e0m0 * hq / 2

F = ufloat(0.0000866, 0.0000001)  # m^2

l = ufloat(0.135, 0.001)  # m
n = 250
R = ufloat(0.7, 0.1)  # ohm
Us = ufloat(0.8, 0.1)

# nu, U = np.genfromtxt('Messung1.txt', unpack=True)
#
# x = np.linspace(0, 100, 2)
# y = [2.3, 2.3]
# y = np.sqrt(y)
# errX = 0.01
# errY = 0.0005
# plt.errorbar(nu, U, xerr=errX, yerr=errY, fmt='none', ecolor='r', label="Messwerte")
# plt.plot(nu, U, 'r.')
# plt.plot(x, y, 'r-', label='$v_0/\sqrt{2}$')
# plt.xlabel('Frequenz $v$/kHz')
# plt.ylabel('Spannung $U$/V')
# # plt.xscale('log')
# # plt.yscale('log')
# plt.xlim(29.8, 40.2)
# # plt.ylim(,)
# plt.tight_layout()
# plt.legend(loc='best')
# # plt.show()
# plt.savefig('Plot1.pdf')

# 35.00  1.55
# 35.20  2.30
# 35.40  1.50

# print('v m p 0')
# print(vm)
# print(vp)
# print(v0)
# Q = v0 / (vp - vm)
# print('Q')
# print(Q)
Mp1 = ufloat(0.01408, 0.0001)
L1 = ufloat(0.156, 0.001)
L1 = ufloat(0.16, 0.001)
Mp2 = ufloat(0.0185, 0.0001)
L2 = ufloat(0.156, 0.001)
Mp3 = ufloat(0.0090, 0.0001)
L3 = ufloat(0.156, 0.001)

rohw1 = ufloat(7400, 10)  # kg/m
rohw2 = ufloat(7800, 100)  # kg/m
rohw3 = ufloat(7240, 10)  # kg/m

Qr1 = Mp1 / (L1 * rohw1)
Qr2 = Mp2 / (L2 * rohw2)
Qr3 = Mp3 / (L3 * rohw3)
# Qr1=11.21
# Qr2=15.77
# Qr3=7.51
N1 = rohw1 / M1 * 2
# N1=2.3*10**28
N2 = rohw2 / M2 * 2
N3 = rohw3 / M3 * 2
print('N')
print(N1)
print(N2)
print(N3)
print('Qreal')
print(Qr1)
print(Qr2)
print(Qr3)

#-------------------------------------------------------------------------------

xib11 = 2 * dR11 * F / (R3 * Qr1)
xib12 = 2 * dR12 * F / (R3 * Qr1)
xib13 = 2 * dR13 * F / (R3 * Qr1)

xib21 = 2 * dR21 * F / (R3 * Qr2)
xib22 = 2 * dR22 * F / (R3 * Qr2)
xib23 = 2 * dR23 * F / (R3 * Qr2)

xib31 = 2 * dR31 * F / (R3 * Qr3)
xib32 = 2 * dR32 * F / (R3 * Qr3)
xib33 = 2 * dR33 * F / (R3 * Qr3)
xib34 = 2 * dR34 * F / (R3 * Qr3)

xib1g = (xib11 + xib12 + xib13) / 3
xib2g = (xib21 + xib22 + xib23) / 3
xib3g = (xib31 + xib32 + xib33 + xib34) / 4

#-------------------------------------------------------------------------------
xia11a = (Ua11 - UBr11) / Us * 4 * l / (omega * muh0 * n**2 * Qr1) * unp.sqrt(R**2 + omega**2 * (muh0 * n**2 / l * F)**2)
xia12a = (Ua12 - UBr12) / Us * 4 * l / (omega * muh0 * n**2 * Qr1) * unp.sqrt(R**2 + omega**2 * (muh0 * n**2 / l * F)**2)
xia13a = (Ua13 - UBr13) / Us * 4 * l / (omega * muh0 * n**2 * Qr1) * unp.sqrt(R**2 + omega**2 * (muh0 * n**2 / l * F)**2)
xia1ga = (xia11a + xia12a + xia13a) / 3
xia21a = (Ua21 - UBr21) / Us * 4 * l / (omega * muh0 * n**2 * Qr2) * unp.sqrt(R**2 + omega**2 * (muh0 * n**2 / l * F)**2)
xia22a = (Ua22 - UBr22) / Us * 4 * l / (omega * muh0 * n**2 * Qr2) * unp.sqrt(R**2 + omega**2 * (muh0 * n**2 / l * F)**2)
xia23a = (Ua23 - UBr23) / Us * 4 * l / (omega * muh0 * n**2 * Qr2) * unp.sqrt(R**2 + omega**2 * (muh0 * n**2 / l * F)**2)
xia2ga = (xia21a + xia22a + xia23a) / 3
xia31a = (Ua31 - UBr31) / Us * 4 * l / (omega * muh0 * n**2 * Qr3) * unp.sqrt(R**2 + omega**2 * (muh0 * n**2 / l * F)**2)
xia32a = (Ua32 - UBr32) / Us * 4 * l / (omega * muh0 * n**2 * Qr3) * unp.sqrt(R**2 + omega**2 * (muh0 * n**2 / l * F)**2)
xia33a = (Ua33 - UBr33) / Us * 4 * l / (omega * muh0 * n**2 * Qr3) * unp.sqrt(R**2 + omega**2 * (muh0 * n**2 / l * F)**2)
xia34a = (Ua34 - UBr34) / Us * 4 * l / (omega * muh0 * n**2 * Qr3) * unp.sqrt(R**2 + omega**2 * (muh0 * n**2 / l * F)**2)
xia3ga = (xia31a + xia32a + xia33a + xia34a) / 4

# selektivverstärker faktor 100
print('test')
print(4 * F * UBr11 / (Qr1 * Us))
#-------------------------------------------------------------------------------
print('xib1')
print(xib11)
print(xib12)
print(xib13)
print(xib1g)
print('xib2')
print(xib21)
print(xib22)
print(xib23)
print(xib2g)
print('xib3')
print(xib31)
print(xib32)
print(xib33)
print(xib34)
print(xib3g)
print('xia1')
print(-xia11a)
print(-xia12a)
print(-xia13a)
print(-xia1ga)
print('xia2')
print(-xia21a)
print(-xia22a)
print(-xia23a)
print(-xia2ga)
print('xia3')
print(-xia31a)
print(-xia32a)
print(-xia33a)
print(-xia34a)
print(-xia3ga)

gj1 = (3 * j1 * (j1 + 1) + (s1 * (s1 + 1) - l1 * (l1 + 1))) / (2 * j1 * (j1 + 1))
gj2 = (3 * j2 * (j2 + 1) + (s2 * (s2 + 1) - l2 * (l2 + 1))) / (2 * j2 * (j2 + 1))
gj3 = (3 * j3 * (j3 + 1) + (s3 * (s3 + 1) - l3 * (l3 + 1))) / (2 * j3 * (j3 + 1))
print('gj')
print(gj1)
print(gj2)
print(gj3)

xic1 = muh0 * muhB**2 * gj1**2 * N1 * j1 * (j1 + 1) / (3 * k * T)
xic2 = muh0 * muhB**2 * gj2**2 * N2 * j2 * (j2 + 1) / (3 * k * T)
xic3 = muh0 * muhB**2 * gj3**2 * N3 * j3 * (j3 + 1) / (3 * k * T)
print('xic')
print(xic1)
print(xic2)
print(xic3)
print('relfehler')
relf = (xic1 - xia1ga) / xic1
print(relf)
# print((xic1-xia1ga)/xic1)
# print((xic2-xia2ga)/xic2)
# print((xic3-xia3ga)/xic3)
# print((xic1-xib1g)/xic1)
# print((xic2-xib2g)/xic2)
# print((xic3-xib3g)/xic3)
# v m p 0
# 35000.0+/-10.0
# 35400.0+/-10.0
# 35200.0+/-10.0
# Q
# 88.0+/-3.1
# Qreal
# 0.1247+/-0.0012
# 0.1603+/-0.0014
# 0.0740+/-0.0013
# xib1
# 0.0002534+/-0.0000028
# 0.0002631+/-0.0000029
# 0.0002436+/-0.0000028
# 0.0002534+/-0.0000026
# xib2
# 0.000419+/-0.000004
# 0.000423+/-0.000004
# 0.000427+/-0.000004
# 0.000423+/-0.000004
# xib3
#(5.40+/-0.25)e-05
#(8.92+/-0.29)e-05
#(8.92+/-0.29)e-05
#(8.21+/-0.28)e-05
#(7.86+/-0.19)e-05
# xia1
# 0.0412+/-0.0018
# 0.0418+/-0.0018
# 0.0379+/-0.0016
# 0.0403+/-0.0011
# xia2
# 0.0813+/-0.0033
# 0.0742+/-0.0028
# 0.0773+/-0.0030
# 0.0776+/-0.0019
# xia3
# 0.00619+/-0.00026
# 0.0086+/-0.0004
# 0.00875+/-0.00034
# 0.0089+/-0.0004
# 0.00810+/-0.00018

# http://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=avogadro+constant
# nist.gov/cgi-bin/cuu/Value?mec2|search_for=electron+mass
# nist.gov/cgi-bin/cuu/Value?k|search_for=boltzman+constant
# http://physics.nist.gov/cgi-bin/cuu/Value?h
# http://physics.nist.gov/cgi-bin/cuu/Value?ukg
# http://www.lenntech.com/calculators/molecular/molecular-weight-calculator.htm
