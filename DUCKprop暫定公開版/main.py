# coding:utf-8

# 改修5.桁に依存しない翼弦長分布のために少し書き換え

import datetime

import os

import numpy
from scipy.optimize import minimize
from scipy import interpolate
import matplotlib.pyplot as plt
import csv
import shutil

from constants_and_settings import *
from represent_with_gamma import *

gamma = numpy.sin(numpy.linspace(0,1,N) * numpy.pi) + 0.01

print('start optimizing gamma...')

power_condition = lambda gamma: -F_P(gamma)+Power if print('power :', F_P(gamma)-Power) or True else numpy.inf
thrust_condition = lambda gamma: F_T(gamma) - Thrust if print('thrust :',F_T(gamma) - Thrust ) or True else numpy.inf
con = (
    #{'type':'eq', 'fun':power_condition},
    {'type':'eq', 'fun':thrust_condition}
)
result = minimize(
    # minusT,
    F_P,
    x0=gamma,
    constraints=con,
    options={"maxiter":1e8,"ftol":0.01},
    method='SLSQP'
)
print('optimization completed')
gamma = result.x
print(result)

CL = F_CL(gamma)
CD = F_CD(gamma)
V  = F_V(gamma)
dL = F_L(gamma)
dD = F_D(gamma)
Re  = F_V(gamma) * F_C(gamma) / NU
phi = numpy.rad2deg(F_phi(gamma)) # 誘導迎え角
alpha = F_alpha(gamma) # 有効迎え角
theta = phi + alpha # 取付迎角
thrust = F_T(gamma)
power=F_P(gamma)
eta = thrust * U_inf / power
chord = F_C(gamma)

print(f'thrust: {thrust}')
print(f'power : {F_P(gamma)}')
print(f'eta   : {eta}')

table = numpy.concatenate(
    [r, chord, theta, phi, alpha, gamma, Re, CL, CD, V, dL, dD]
).reshape((12,-1 )).T
header = '\n'.join(
    [
        f'RHO   : {RHO}',
        f'NU    : {NU}',
        f'B     : {B}',
        f'U_inf : {U_inf}',
        f'RPM   : {RPM}',
        f'power : {power}',
        f'R_out : {R_out}',
        f'R_in  : {R_in}',
        f'N     : {N}',
        f'M     : {M}',
        f'dt    : {dt}',
        f'thrust: {thrust}',
        f'eta   : {eta}',
        f'optsucess:{result.message}'
    ]
) + '\n' + ','.join(
    [
        '半径r[m]',
        '翼弦長chord[m]',
        '取付角theta[deg]',
        '流入角phi[deg]',
        '迎角alpha[deg]',
        '循環gamma',
        'Re',
        '揚力係数CL',
        '抗力係数CD',
        '流入速度V',
        '局所揚力dL',
        '局所抗力dD'
    ]
)

now = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')

folder="prop_result_"+now
os.makedirs(folder)

filename = folder+'/'+'result_at_'+now+'.txt'
numpy.savetxt(filename, table, delimiter=',', header=header, comments='')

r_spline = numpy.arange(0, R_out, 0.001)
chord_spline = interpolate.interp1d(
    r, chord, kind='linear', fill_value='extrapolate'
)
chord_spline_list=chord_spline(r_spline)

theta_spline = interpolate.interp1d(
    r, theta, kind='linear', fill_value='extrapolate'
)
theta_spline_list=theta_spline(r_spline)

header = ['半径r[mm]','翼弦長chord[mm]','取付角theta[deg]']
filename = folder+'/'+'for_design_at'+now+'.csv'
with open(filename,'w',newline="") as f:
    writer=csv.writer(f)
    writer.writerow(header)
    for i in range(len(r_spline)):
        writer.writerow(
            [r_spline[i], chord_spline_list[i], theta_spline_list[i]]
            )

fig=plt.figure(figsize=(12,8.4))

plt.subplot(3,3,1)
plt.plot(r, chord)
plt.title('Chord(chord length)')

plt.subplot(3,3,2)
plt.plot(r, theta)
plt.title('Theta(mounting anle)')

plt.subplot(3,3,3)
plt.plot(r, gamma)
plt.title('gam(circulation)')

plt.subplot(3,3,4)
plt.plot(r, Re)
plt.title('Re')

plt.subplot(3,3,5)
plt.plot(r, CL)
plt.title('Cl')

plt.subplot(3,3,6)
plt.plot(r, CD)
plt.title('Cd')

plt.subplot(3,3,7)
plt.plot(r, V)
plt.title('V')

plt.subplot(3,3,8)
plt.plot(r, dL)
plt.title('dL')

plt.subplot(3,3,9)
plt.plot(r, dD)
plt.title('dD')

fig.savefig(folder+'/'+now+'_propfig.png')

fig=plt.figure(figsize=(10,13))
fig.subplots_adjust(hspace=0.2, wspace=0.3)

plt.subplot(4,2,1)
plt.plot(r, chord)
plt.ylabel('Chord(chord length)[m]')

plt.subplot(4,2,2)
plt.plot(r, theta)
plt.ylabel('Theta(mounting anle)[deg]')

plt.subplot(4,2,3)
plt.plot(r, gamma)
plt.ylabel('gam(circulation)[m^2/s]')

plt.subplot(4,2,4)
plt.plot(r, Re)
plt.ylabel('Re(Reynolds number)')

plt.subplot(4,2,5)
plt.plot(r, CL)
plt.ylabel('Cl(lift coefficient)')

plt.subplot(4,2,6)
plt.plot(r, CD)
plt.ylabel('Cd(drag coefficient)')

plt.subplot(4,2,7)
plt.plot(r, V)
plt.ylabel('V(flow velocity)[m/s]')
plt.xlabel('r(distance from center of rotation)[m]')

plt.subplot(4,2,8)
plt.plot(r, alpha)
plt.ylabel('alpha(angle of attack)[deg]')
plt.xlabel('r(distance from center of rotation)[m]')


fig.savefig(folder+'/'+now+'_propfig_thesis.png')