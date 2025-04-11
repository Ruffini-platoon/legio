# coding: utf-8

#改修1.ペラ根の後付けでの改変を前提としてペラ桁に関する部分を削除

import numpy
from scipy import interpolate

def sep_blade(R_out, R_in, N):
    delta_B = (R_out - R_in) / N # 翼素の長さ[m]
    r = numpy.linspace(R_in, R_out-delta_B, N) + delta_B/2 # 各翼素の中心の回転半径
    # 各翼素の中心の位置ベクトル
    CP = numpy.array([
        numpy.array([0, r[i], 0]) for i in range(N)
    ])
    return CP, delta_B, r

def sep_vortex(CP, delta_B, r, U_inf, Omega, B, N, M, dt):
    r = r - delta_B/2
    r = numpy.concatenate([r, [r[-1]+delta_B]]) # 内外の端が必要なので

    DP = numpy.empty((B, N+1, M+1, 3))
    for bl in range(B):
        for j in range(N+1): # 内外の端を使うため
            for k in range(M+1): # 前端と後端を使うため
                t = k * dt
                theta = 2*numpy.pi/B * bl + Omega * t

                x = - U_inf * t
                y = r[j] * numpy.cos(-theta)
                z = r[j] * numpy.sin(-theta)
                DP[bl,j,k,:] = numpy.array([x, y, z])
    return DP

# 物理定数
RHO = 1.225
NU  = 1.463e-5

# 設計依存の定数
B     = 2 # ブレード数
U_inf = 7.9 # 機速[m/s]
RPM   = 140 # ペラ回転数[rpm]
#Power =  # 出力[W]
Thrust = 20 # 推力[N]
Omega = numpy.deg2rad(RPM*360) / 60 # 回転数[rad/sec]
R_out  = 1.35 # ペラ端回転半径[m]
R_in   = 0.059 # ペラ根回転半径[m]


# 解析のための定数
N  = 50 # ペラ分割数
M  = 50 # 馬蹄渦分割数
dt = 0.05 # 馬蹄渦分割時間幅

# 各代表点の座標
CP, delta_B, r = sep_blade(R_out, R_in, N)
DP = sep_vortex(CP, delta_B, r, U_inf, Omega, B, N, M, dt)
