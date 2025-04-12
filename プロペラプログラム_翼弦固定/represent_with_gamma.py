# coding:utf-8

# 3.2-4節の内容

import numpy

from constants_and_settings import *
# from minimize_cd import f_CD, f_CL, f_C, f_alpha, f_CLRe
from minimize_cd import dae51

# RHO, NU は物理定数
# N, M, B, CP, CD, delta_B, r, U_inf, Omega は定数

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

# 対気速度のx成分(式6)
F_Up = lambda gamma: U_inf - u_ind(gamma)

# 対気速度のz成分(式5)
F_Ut = lambda gamma: Omega * r - w_ind(gamma)

# 対気速度の大きさ(式7, 18)
F_V  = lambda gamma: numpy.sqrt(F_Up(gamma)**2+F_Ut(gamma)**2)

# 誘導迎え角(式8, 19)
F_phi = lambda gamma: numpy.arctan(F_Up(gamma)/F_Ut(gamma))

# 局所揚力(式9, 20)
F_L = lambda gamma: RHO * F_V(gamma) * delta_B * gamma

chord = (dia_beam+RIB_MARGIN)/Thickness
F_Re  = lambda gamma: chord*F_V(gamma)/NU
F_CL  = lambda gamma: F_L(gamma) / (RHO*F_V(gamma)**2*chord*delta_B/2)
F_alpha = lambda gamma: dae51.alpha_ReCL(F_Re(gamma), F_CL(gamma))
F_CD  = lambda gamma: dae51.CD_ReCL(F_Re(gamma), F_CL(gamma))

# gammaのみの関数としてのdD(式10, 25)
F_D = lambda gamma: (1/2) * RHO * F_V(gamma)**2 * chord * delta_B * F_CD(gamma)

# gamma のみの関数としての推力(式12, 27)
def F_T(gamma):
    phi = F_phi(gamma)
    T = B * numpy.sum(
        F_L(gamma) * numpy.cos(phi) - F_D(gamma) * numpy.sin(phi)
    )

    return T

minusT = lambda gamma: -F_T(gamma)

# gammaのみの関数としてのパワー(式13, 28)
def F_P(gamma):
    phi = F_phi(gamma)
    P = B * Omega * numpy.sum(
        ( F_D(gamma) * numpy.cos(phi) + F_L(gamma) * numpy.sin(phi) ) * r
    )
    return P
