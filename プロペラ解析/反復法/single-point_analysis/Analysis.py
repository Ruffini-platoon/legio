# coding: utf-8

import numpy
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import minimize

from Load_Prop import large_table , delta_B , CP , N
from load_and_interpolate import dae51

r = large_table[:,0]
chord = large_table[:,1]
theta = large_table[:,2]

#解析したい点を入力
U_inf = 9 # 機速[m/s]
RPM   = 140 # ペラ回転数[rpm]

# 物理定数
RHO = 1.225
NU  = 1.463e-5

# 設計依存の定数
B     = 2 # ブレード数
Omega = numpy.deg2rad(RPM*360) / 60 # 回転数[rad/sec]

# 解析のための定数
#ペラ分割数NはLoad_Propに
M  = 50 # 馬蹄渦分割数
dt = 0.05 # 馬蹄渦分割時間幅

def sep_vortex(CP, delta_B, r, U_inf, Omega, B, N, M, dt):
    r = r - delta_B/2
    r = numpy.concatenate([r, [r[-1]+delta_B]]) # 内外の端が必要なので

    DP = numpy.empty((B, N+1, M+1, 3))
    for bl in range(B):
        for j in range(N+1): # 内外の端を使うため
            for k in range(M+1): # 前端と後端を使うため
                t = k * dt
                thita = 2*numpy.pi/B * bl + Omega * t

                x = - U_inf * t
                y = r[j] * numpy.cos(-thita)
                z = r[j] * numpy.sin(-thita)
                DP[bl,j,k,:] = numpy.array([x, y, z])
    return DP

# 各代表点の座標
DP = sep_vortex(CP, delta_B, r, U_inf, Omega, B, N, M, dt)

# a, bを結ぶ線分上の単位強度の循環が点Oに引き起こす流速(式1)
def f(a, b):
    # a, bは一次元配列
    l = b - a
    al = numpy.cross(a, l)  # aとlの外積
    e_a = a / numpy.linalg.norm(a)
    e_b = b / numpy.linalg.norm(b)
    if(numpy.linalg.norm(al)==0):
        f=0
    else:
        f = 1/(4*numpy.pi) * (al/numpy.linalg.norm(al)**2) * numpy.dot((e_b - e_a), l)
    return f

# j番目の馬蹄渦(単位強度とする)がi番目の翼素に引き起こす流速(式2)
def ind_i_j(CP, DP):
    # CP[i](i番目の翼素上の制御点の位置ベクトル)は関数外で定義
    # DP[bl,j,k](bl番目のブレードからのj番目の放出渦のk番目の分割点)貼関数外で定義
    # N(ブレードの分割数)、B(ブレードの数-1)、M(馬蹄渦の分割数-1)は関数外で定義されているものとする
    print()
    print('start computing biot-savart coefficient matrix')
    a = numpy.array(
        [
            [
                [
                    [
                        DP[bl, j, k] - CP[i]
                        for k in range(M)
                    ]
                    for j in range(N)
                ]
                for i in range(N)
            ]
            for bl in range(B)
        ]
    )
    b = numpy.array(
        [
            [
                [
                    [
                        DP[bl, j, k+1] - CP[i]
                        for k in range(M)
                    ]
                    for j in range(N)
                ]
                for i in range(N)
            ]
            for bl in range(B)
        ]
    )
    a_dash = numpy.array(
        [
            [
                [
                    [
                        DP[bl, j+1, k] - CP[i]
                        for k in range(M)
                    ]
                    for j in range(N)
                ]
                for i in range(N)
            ]
            for bl in range(B)
        ]
    )
    b_dash = numpy.array(
        [
            [
                [
                    [
                        DP[bl, j+1, k+1] - CP[i]
                        for k in range(M)
                    ]
                    for j in range(N)
                ]
                for i in range(N)
            ]
            for bl in range(B)
        ]
    )
    print('  complete setting relative position')
    X = numpy.empty((N,N))
    Y = numpy.empty((N,N))
    Z = numpy.empty((N,N))

    print()
    print('  computing coefficient matrix')
    for i in range(N):
        print(f'    {100*i/N}%...')
        for j in range(N):
            s1 = numpy.sum(
                [
                    numpy.sum(
                        [
                            f(a[bl,i,j,k],b[bl,i,j,k])
                            for k in range(M)
                        ],
                        axis=0
                    )
                    for bl in range(B)
                ],
                axis=0
            )
            s2 = numpy.sum(
                [
                    f(a[bl,i,j,1], a_dash[bl,i,j,1])
                    for bl in range(B)
                ],
                axis=0
            )
            s3 = numpy.sum(
                [
                    numpy.sum(
                        [
                            f(a_dash[bl,i,j,k],b_dash[bl,i,j,k])
                            for k in range(M)
                        ],
                        axis=0
                    )
                    for bl in range(B)
                ],
                axis=0
            )
            X[i,j], Y[i,j], Z[i,j] = -s1 + s2 + s3

    print('  complete')

    return X, Y, Z

X, Y, Z = ind_i_j(CP, DP)
# 以下X, Zも定数として与えられているとする

# ブレード全体の循環が引き起こすx方向の誘導速度(式3)
def u_ind(gamma):
    # X[i,j]はj番目の馬蹄渦(単位強度)がi番目の翼素に引き起こす流速のX成分
    # gamma[j]はj番目の馬蹄渦の(実際の)循環
    u = numpy.array(
        [
            numpy.sum(
                [
                    X[i,j] * gamma[j]
                    for j in range(N)
                ]
            )
            for i in range(N)
        ]
    )
    return u

# ブレード全体の循環が引き起こすz方向の誘導速度(式4)
def w_ind(gamma):
    # Z[i,j]はj番目の馬蹄渦(単位強度)がi番目の翼素に引き起こす流速のZ成分
    # gamma[j]はj番目の馬蹄渦の(実際の)循環
    w = numpy.array(
        [
            numpy.sum(
                [
                    Z[i,j] * gamma[j]
                    for j in range(N)
                ]
            )
            for i in range(N)
        ]
    )
    return w




#推力と効率をgammaの関数にする

# 対気速度のx成分
F_Up = lambda gamma: U_inf - u_ind(gamma)

# 対気速度のz成分
F_Ut = lambda gamma: Omega * r - w_ind(gamma)

# 対気速度の大きさ
F_V  = lambda gamma: numpy.sqrt(F_Up(gamma)**2+F_Ut(gamma)**2)

# 流入角度
F_phi = lambda gamma: numpy.arctan(F_Up(gamma)/F_Ut(gamma))

# 有効迎角
F_alpha = lambda gamma: theta - numpy.rad2deg(F_phi(gamma))

# Re数
F_Re = lambda gamma: F_V(gamma)*chord / NU

#揚力係数
#F_CL = lambda gamma: dae51.CL_alphaRe(F_alpha(gamma),F_Re(gamma)) #使用する翼型が二つ以下ならlanmbdaで定義しても良い。
    
#半径位置によって参照する翼型データを場合分け
def F_CL( gamma ):
    for gam, rad in zip(gamma, r):
        if rad < 1.45:
            return dae51.CL_alphaRe(F_alpha(gamma),F_Re(gamma))
        else:
            return dae51.CL_alphaRe(F_alpha(gamma),F_Re(gamma)) 

#抗力係数
#F_CD = lambda gamma: dae51.CD_alphaRe(F_alpha(gamma),F_Re(gamma)) #使用する翼型が二つ以下ならlanmbdaで定義しても良い。
    
#半径位置によって参照する翼型データを場合分け
def F_CD( gamma ):
    for gam, rad in zip(gamma, r):
        if rad < 1.45:
            return dae51.CD_alphaRe(F_alpha(gamma),F_Re(gamma))
        else:
            return dae51.CD_alphaRe(F_alpha(gamma),F_Re(gamma))

#局所揚力
F_L = lambda gamma: (1/2)*delta_B*RHO*chord*F_CL(gamma)*(F_V(gamma))**2

#局所抗力
F_D = lambda gamma: (1/2)*delta_B*RHO*chord*F_CD(gamma)*(F_V(gamma))**2

#循環
F_gamma = lambda gamma: F_L(gamma) / (RHO * F_V(gamma) * delta_B)

# gamma のみの関数としての推力
def F_T(gamma):
    phi = F_phi(gamma)
    T = B * numpy.sum(
        F_L(gamma) * numpy.cos(phi) - F_D(gamma) * numpy.sin(phi)
    )

    return T

# gammaのみの関数としてのぱわーーーー
def F_P(gamma):
    phi = F_phi(gamma)
    P = B * Omega * numpy.sum(
        ( F_D(gamma) * numpy.cos(phi) + F_L(gamma) * numpy.sin(phi) ) * r
    )
    return P

#効率
F_eta = lambda gamma: U_inf * F_T(gamma) / F_P(gamma)



'''
#初期値（誘導速度を考慮しない揚力分布より求まる循環分布）を求める。

#初期大気速度
fV = numpy.sqrt(U_inf**2 + (Omega * r)**2)

#初期流入角度
fphi = numpy.arctan(U_inf / (Omega * r))

#初期有効迎角
falpha = theta - numpy.rad2deg(fphi)

#初期Re数
fRe = fV * chord / NU

#初期揚力係数
fCL = dae51.CL_alphaRe(falpha , fRe)

#初期抗力係数
fCD = dae51.CD_alphaRe(falpha , fRe)

#初期局所揚力
fL = (1/2)*delta_B*RHO*chord*fCL*(fV)**2

#初期局所抗力
fD = (1/2)*delta_B*RHO*chord*fCD*(fV)**2

#初期推力
fT = B * numpy.sum(fL * numpy.cos(fphi) - fD * numpy.sin(fphi))

#初期ぱわーーーーーーー
fP = B * Omega * numpy.sum(( fD * numpy.cos(fphi) + fL * numpy.sin(fphi) ) * r)

#初期効率
feta = U_inf * fT / fP

#初期循環
fgamma = fL / (RHO * fV * delta_B)


#当初試みていた反復法の残骸。途中でminimizeの方が良いことに気づいてやめちゃった。一応残しておく。

#反復
gamma = fgamma #初期値
T = fT #初期値
T_pre = 1 #仮置き
a = 1


print('start iteration')
while numpy.abs( T - T_pre ) > 1 :
    T_pre = T #T_preを更新
    T = F_T(gamma) #Tを更新
    eta = F_eta(gamma)
    print('Iter' , a , ': ' , T , 'N, ', eta , '%')
    gamma = F_gamma(gamma) #gammaを更新
    a = a + 1
print('complete')
print('解析結果 @ U_inf = ' ,U_inf,' RPM = ',RPM,)
print(F_T(gamma),' N')
print(F_eta(gamma),' %')
'''


#パソコンさんにgammaとF_gammaが(ほぼ)等しくなるgammaを探してもらう
print()
print('start searching gamma')

#目的関数
Object = lambda gamma: numpy.linalg.norm( F_gamma(gamma) - (gamma) ) if print('thrust :',F_T(gamma) ) or True else numpy.inf

#初期値
ffgamma = numpy.sin(numpy.linspace(0,1,N) * numpy.pi) + 0.01
#ffgamma = fgamma #やめたほうがいい

if numpy.isnan(F_T(ffgamma)):#初期値で推力がnanになったらあきらめる
    print('解析結果 @ U_inf = ' ,U_inf,' RPM = ',RPM,)
    print('Analysis failed')

else:
    result = minimize(
            Object,
            x0=ffgamma,
            options={"maxiter":1e8,"ftol":1e-2},#ftolは1e-3くらいが良い
            method='SLSQP'
        )
    print('complete')
    print()

    #結果の表示
    print('解析結果 @ U_inf = ' ,U_inf,' RPM = ',RPM,)
    print(F_T(result.x),' N')
    print(F_eta(result.x),' %')
