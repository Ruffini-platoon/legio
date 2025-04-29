# coding:utf-8

import os
import re
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy

from Load_Prop import large_table , delta_B , CP , N
from load_and_interpolate import dae51
from Analysis import *

#メッシュグリッドの作成
U_inf_vals = numpy.linspace(4 , 10 , 4) #適宜変更
RPM_vals = numpy.linspace(80 , 170 , 4) #適宜変更(80~170にする)

U_grid, RPM_grid = numpy.meshgrid(U_inf_vals, RPM_vals)

result_T = numpy.zeros_like(U_grid)
result_eta = numpy.zeros_like(U_grid)

for i in range(U_grid.shape[0]):
    for j in range(U_grid.shape[1]):
        print()
        print('start analysing @ U_inf = ' ,U_grid[i,j],' RPM = ',RPM_grid[i,j],)
        result_T[i,j] ,result_eta[i,j] = Robust_analysis_at( U_grid[i, j], RPM_grid[i, j] )



#推力分布をグラフに
fig = plt.figure(figsize=(8, 6))
contour = plt.contourf(U_grid, RPM_grid, result_T, levels=50, cmap='viridis')
plt.colorbar(contour)

plt.xlabel('U_inf')
plt.ylabel('RPM')
plt.title('Contour Plot of Thrust')
fig.savefig('Thrust robustness analysis.png')
plt.show()


#効率の変動をグラフに
fig = plt.figure(figsize=(8, 6))
contour = plt.contourf(U_grid, RPM_grid, result_eta, levels=50, cmap='viridis')
plt.colorbar(contour)

plt.xlabel('U_inf')
plt.ylabel('RPM')
plt.title('Contour Plot of efficiency')
fig.savefig('efficiency robustness analysis.png')
plt.show()