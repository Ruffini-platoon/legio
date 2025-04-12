# coding: utf-8

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

def interpolate_dia_beam(dia_beam_root, dia_beam_tip, r_beam_root, r_beam_tip, r):
    fun = interpolate.interp1d(
        [r_beam_root, r_beam_tip],
        [dia_beam_root, dia_beam_tip],
        kind='linear', bounds_error=False, fill_value='extrapolate'
    )
    return fun(r)

# 物理定数
RHO = 1.225
NU  = 1.463e-5

# 設計依存の定数
B     = 2 # ブレード数
U_inf = 7.9 # 機速[m/s]
RPM   = 140 # ペラ回転数[rpm]
Power = 216 # 出力[W]
Thrust =  23 # 推力[N]
Omega = numpy.deg2rad(RPM*360) / 60 # 回転数[rad/sec]
R_out  = 1.45 # ペラ端回転半径[m]
R_in   = 0.06 # ペラ根回転半径[m]
dia_beam_root = 0.019748 # 桁根の外径
dia_beam_tip=0.00636 # 桁端の外径
r_beam_root = 0.006 # 桁根の回転半径
beam_length = 1.25 # 桁長さ
r_beam_tip = r_beam_root + beam_length
Thickness = 0.0938 # 最大翼厚(/chord)
RIB_MARGIN = 0.006 # リブの翼厚の余白


# 解析のための定数
N  = 50 # ペラ分割数
M  = 50 # 馬蹄渦分割数
dt = 0.05 # 馬蹄渦分割時間幅

# 各代表点の座標
CP, delta_B, r = sep_blade(R_out, R_in, N)
DP = sep_vortex(CP, delta_B, r, U_inf, Omega, B, N, M, dt)
dia_beam = interpolate_dia_beam(
    dia_beam_root,
    dia_beam_tip,
    r_beam_root,
    r_beam_tip,
    r
)
