# coding:utf-8

# 3.5, 3.6節の内容の実装
# 規格化はしていない

import os
import re

import matplotlib.pyplot as plt

import numpy
from scipy import interpolate
from scipy import optimize
import matplotlib.pyplot as plt


from constants_and_settings import *

numpy.set_printoptions(
    threshold = numpy.inf,
    linewidth = 1e3
)

# 翼型の解析データを置くディレクトリ
AIRFOIL_DIR = './airfoil'


class AirFoil:
    def __init__(self, airfoil_name):
        # データの読み込み
        large_table = self.load_files(airfoil_name)
        Re = large_table[:,0]
        alpha = large_table[:,1]
        CL = large_table[:,2]
        CD = large_table[:,3]

        Re_grid = (60_000, 300_000, 50)
        CL_grid = (0.2, 0.9, 50)
        Re = numpy.linspace(*Re_grid)
        CL = numpy.linspace(*CL_grid)

        x, y = numpy.meshgrid(Re,CL)

        alphaReCL_grid = interpolate.griddata(
            (large_table[:,0], large_table[:,2]),
            large_table[:,1],
            (x,y),
            method='linear'
        )
        alphaReCL=interpolate.interp2d(
                Re, CL, alphaReCL_grid,
                kind="linear",
                bounds_error = False
            )
        self.alpha_ReCL= lambda Re,CL:numpy.diag(alphaReCL(Re,CL))

        #self.alpha_ReCL = self.refine_interp2d(
        #    interpolate.interp2d(
        #        Re, CL, alphaReCL_grid,
        #        kind="linear",
        #        bounds_error = False
        #    )
        #)

        CDReCL_grid = interpolate.griddata(
            (large_table[:,0], large_table[:,2]),
            large_table[:,3],
            (x,y),
            method='linear'
        )
        CDReCL=interpolate.interp2d(
                Re, CL, CDReCL_grid,
                kind="linear",
                bounds_error = False
            )
        self.CD_ReCL= lambda Re,CL:numpy.diag(CDReCL(Re,CL))

        #self.CD_ReCL = self.refine_interp2d(
        #    interpolate.interp2d(
        #        Re, CL, CDReCL_grid,
        #        kind="linear",
        #        bounds_error = False
        #    )
        #)

        print('airfoil data successfully loaded')
        return

    @staticmethod
    def load_files(airfoil_name):
        print()
        # ファイルの検索
        files   = os.listdir('./'+AIRFOIL_DIR+'/'+airfoil_name)
        for filename in files:
            if airfoil_name+'_T1' not in filename:
                print('the data of different airfoil included')
                print('check the contents of the folder and try again')
                exit()
        files = ['./'+AIRFOIL_DIR+'/'+airfoil_name+'/'+filename for filename in files]
        print(f'start loading {airfoil_name}...')

        large_table = None

        fig=plt.figure()
        # データの読み込み
        for row_n, filename in enumerate(files):
            match = re.search('Re\d\.\d{3}', filename)
            Re    = numpy.float64(match.group()[2:]) * 1e6
            print(f'    loading file of Re={Re}')
            table = numpy.loadtxt(filename, skiprows = 10, usecols = [0, 1, 2], delimiter=',')
            # Re, alpha, Cl, Cd の形の配列にする
            Re_col = numpy.repeat(Re, table.shape[0]).reshape((-1, 1))
            table = numpy.concatenate([Re_col, table], axis=1)

            plt.plot(table[:,3],table[:,2])

            # large_table 縦に結合する
            if type(large_table) == type(None):
                large_table = table
            else:
                large_table = numpy.concatenate([large_table, table], axis=0)

        fig.savefig('polar.png')

        print('complete')
        print()
        return large_table
    
    """
    @staticmethod
    def refine_interp2d(interp2d):
        def fun(x, y, mesh=False):
            '''
            scipy.interpolate.interp2d の返り値の形を整えた関数

            Parameters
            ----------
            x : float or numpy.ndarray
            y : float or numpy.ndarray
            mesh : boolean
            '''
            try:
                if type(x) != numpy.ndarray and type(y) != numpy.ndarray:
                    # 引数が二つともスカラー -> スカラーで返す(もとは長さ1の配列)
                    return interp2d(x, y)[0]
                elif type(x) != numpy.ndarray and type(y) == numpy.ndarray:
                    # 第2引数のみベクトル -> 一次元配列で返す(もとは列ベクトル(二次元配列))
                    return interp2d(x, y).reshape(-1)
                elif type(x) == numpy.ndarray and type(y) != numpy.ndarray:
                    # 第1引数のみベクトル -> 一次元配列で返す(もとと同じ)
                    return interp2d(x, y)
                elif type(x) == numpy.ndarray and type(y) == numpy.ndarray:
                    # 両引数がベクトル
                    if mesh:
                        # mesh 指定がある -> 両ベクトルが張るメッシュの値を持つ二次元配列を返す(もとのまま)
                        return interp2d(x, y)
                    else:
                        # mesh 指定なし -> 両ベクトルの同じ位置の2数に対する値を並べた一次元配列を返す
                        if len(x) == len(y):
                            # 正方行列の対角成分を取る
                            return numpy.diag(interp2d(x, y))
                        else:
                            raise ValueError('two vectors must have the same sizes')
                else:
                    raise TypeError('two arguments must be two vectors, two scalars, or a vector and a scalar')
            except ValueError as e:
                print(f'x: {x}')
                print(f'y: {y}')
                raise e
        return fun
    """

dae51 = AirFoil('DAE51')