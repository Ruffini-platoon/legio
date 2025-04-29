# coding:utf-8

import os
import re

import numpy
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator, griddata
from scipy import optimize
import matplotlib.pyplot as plt

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

        Re_grid = ( 0 , 3e5 , 50)
        alpha_grid = (-5, 20, 50)
        Re = numpy.linspace(*Re_grid)
        alpha = numpy.linspace(*alpha_grid)

        x, y = numpy.meshgrid(alpha,Re)

        self.CL_alphaRe = lambda alpha, Re: interpolate.griddata(
            (large_table[:, 1] , large_table[:, 0]),
            large_table[:, 2],  
            (alpha, Re),
            method='cubic' #論文ではスプライン補間だったため
            )

        self.CD_alphaRe = lambda alpha, Re: interpolate.griddata(
            (large_table[:,1] , large_table[:,0]), 
            large_table[:,3],
            (alpha,Re),
            method='cubic'
        )
        
        self.Re= Re
        self.alpha= alpha

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
            table = numpy.loadtxt(filename, skiprows = 11, usecols = [0, 1, 2], delimiter=',')
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
    
    

dae51 = AirFoil('DAE51')