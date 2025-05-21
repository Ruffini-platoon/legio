# coding:utf-8

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy

from Load_Prop import large_table , delta_B , CP , N
from load_and_interpolate import dae51
from Analysis import *

#メッシュグリッドの作成
U_inf_vals = numpy.linspace(2 , 9 , 8) #適宜変更
RPM_vals = numpy.linspace(90 , 160 , 8) #適宜変更(80~170にする)。0を含めないこと。

U_grid, RPM_grid = numpy.meshgrid(U_inf_vals, RPM_vals)

result_T = numpy.zeros_like(U_grid)
result_P = numpy.zeros_like(U_grid)
result_eta = numpy.zeros_like(U_grid)
result_norm = numpy.zeros_like(U_grid)

for i in range(U_grid.shape[0]):
    for j in range(U_grid.shape[1]):

        print()
        print('start analysing @ U_inf = ' ,U_grid[i,j],' RPM = ',RPM_grid[i,j],)

        result_T[i,j], result_P[i,j], result_eta[i,j], result_norm[i,j] = Robust_analysis_at( U_grid[i, j], RPM_grid[i, j] )

        if result_T[i,j] < 0: #TやPが負だと効率がおかしなことになるので負の場合は除いておく。(nanはplot時に除かれる)
            result_eta[i,j] = numpy.nan
        else:
            if result_P[i,j] < 0:
                result_eta[i,j] = numpy.nan




#折れ線グラフの重ね合わせ
fig = plt.figure(figsize=(12, 10), constrained_layout=True)
# U_infに対するTの変化
plt.subplot(2,2,1)

for i in range(RPM_vals.shape[0]):
    plt.plot(U_inf_vals, result_T[i, :], label=f'RPM = {RPM_vals[i]:.0f}', linewidth=1, marker='o' ,markersize=3)

plt.xlabel('U_inf [m/s]')
plt.ylabel('Thrust T [N]')
plt.title('Thrust vs U_inf')
plt.legend()
plt.grid(True)
#plt.tight_layout() #一つずつ出力していた頃の名残(以下同)
#fig.savefig('Thrust robustness analysis(2D).png')
#plt.show()

# U_infに対するPの変化
plt.subplot(2,2,2)

for i in range(RPM_vals.shape[0]):
    plt.plot(U_inf_vals, result_P[i, :], label=f'RPM = {RPM_vals[i]:.0f}', linewidth=1, marker='o' ,markersize=3)

plt.xlabel('U_inf [m/s]')
plt.ylabel('Power')
plt.title('Power vs U_inf')
plt.legend()
plt.grid(True)
#plt.tight_layout()
#fig.savefig('Power robustness analysis(2D).png')
#plt.show()

# U_infに対するetaの変化
plt.subplot(2,2,3)

for i in range(RPM_vals.shape[0]):
    plt.plot(U_inf_vals, result_eta[i, :], label=f'RPM = {RPM_vals[i]:.0f}', linewidth=1, marker='o' ,markersize=3)

plt.xlabel('U_inf [m/s]')
plt.ylabel('Efficiency')
plt.title('Efficiency vs U_inf')
plt.legend()
plt.grid(True)
#plt.tight_layout()
#fig.savefig('Efficiency robustness analysis(2D).png')
#plt.show()

# U_infに対するetaの変化
plt.subplot(2,2,4)

for i in range(RPM_vals.shape[0]):
    plt.plot(U_inf_vals, result_norm[i, :], label=f'RPM = {RPM_vals[i]:.0f}', linewidth=1, marker='o' ,markersize=3)

plt.xlabel('U_inf [m/s]')
plt.ylabel('norm')
plt.title('norm vs U_inf')
plt.legend()
plt.grid(True)
#plt.tight_layout()

fig.savefig('Robustness analysis(2D).png')
plt.show()


#等位線プロット(個人的に分かりづらい)
fig = plt.figure(figsize=(12, 10), constrained_layout=True)
#推力分布をグラフに
plt.subplot(2,2,1)
contour = plt.contourf(U_grid, RPM_grid, result_T, levels=50, cmap='viridis')
plt.colorbar(contour)

plt.xlabel('U_inf')
plt.ylabel('RPM')
plt.title('Contour Plot of Thrust')
#fig.savefig('Thrust robustness analysis(3D).png')
#plt.show()

#ぱわーの変動をグラフに
plt.subplot(2,2,2)
contour = plt.contourf(U_grid, RPM_grid, result_P, levels=50, cmap='viridis')
plt.colorbar(contour)

plt.xlabel('U_inf')
plt.ylabel('RPM')
plt.title('Contour Plot of Power')
#fig.savefig('Power robustness analysis(3D).png')
#plt.show()

#効率の変動をグラフに
plt.subplot(2,2,3)
contour = plt.contourf(U_grid, RPM_grid, result_eta, levels=50, cmap='viridis')
plt.colorbar(contour)

plt.xlabel('U_inf')
plt.ylabel('RPM')
plt.title('Contour Plot of efficiency')
#fig.savefig('Efficiency robustness analysis(3D).png')
#plt.show()

#効率の変動をグラフに
plt.subplot(2,2,4)
contour = plt.contourf(U_grid, RPM_grid, result_norm, levels=50, cmap='viridis')
plt.colorbar(contour)

plt.xlabel('U_inf')
plt.ylabel('RPM')
plt.title('Contour Plot of norm')
fig.savefig('Robustness analysis(3D).png')
plt.show()