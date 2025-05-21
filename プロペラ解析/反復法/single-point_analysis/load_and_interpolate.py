# coding:utf-8

import os
import re

import numpy
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator, griddata, RBFInterpolator
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

        Re_grid = ( 10000 , 3e5 , 50)
        alpha_grid = (-10, 20, 50)
        Re = numpy.linspace(*Re_grid)
        alpha = numpy.linspace(*alpha_grid)

        x, y = numpy.meshgrid(alpha,Re)

        '''        
        #griddata補間(信頼できるが外挿できない)
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
        '''

        #CLを補間(外挿あり)
        CL_alphaRe_grid = griddata(
            (large_table[:, 1] , large_table[:, 0]),
            large_table[:,2],
            (x, y),
            method='cubic'
        )

        CL_alphaRe_Reg = RegularGridInterpolator(
            (alpha, Re),
            CL_alphaRe_grid.T,
            bounds_error=False,
            fill_value=None
        )

        self.CL_alphaRe = lambda alpha, Re: CL_alphaRe_Reg( numpy.column_stack([alpha , Re]) )

        #CDを補間(外挿あり)
        CD_alphaRe_grid = griddata(
            (large_table[:, 1] , large_table[:, 0]),
            large_table[:,3],
            (x, y),
            method='cubic'
        )

        CD_alphaRe_Reg = RegularGridInterpolator(
            (alpha, Re),
            CD_alphaRe_grid.T,
            bounds_error=False,
            fill_value=None
        )

        self.CD_alphaRe = lambda alpha, Re: CD_alphaRe_Reg( numpy.column_stack([alpha , Re]) )

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
            plt.title(airfoil_name)

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
GEMINI = AirFoil('GEMINI')
ge51_25 = AirFoil('ge51_25')
ge51_50 = AirFoil('ge51_50')
ge51_75 = AirFoil('ge51_75')


#補間確認用
def plot_contour(airfoil, quantity='CL', levels=20):
    alpha = numpy.linspace(-30 ,40 , 500)
    Re = numpy.linspace(10000 , 300000 , 500)
    x, y = numpy.meshgrid(alpha, Re)

    if quantity == 'CL':
        z = airfoil.CL_alphaRe(x.flatten(), y.flatten())
    elif quantity == 'CD':
        z = airfoil.CD_alphaRe(x.flatten(), y.flatten())
    else:
        raise ValueError("quantity must be 'CL' or 'CD'")

    z = z.reshape(x.shape)

    plt.figure(figsize=(8, 6))
    cp = plt.contourf(x, y, z, levels=levels, cmap='viridis')
    plt.colorbar(cp, label=quantity)
    plt.xlabel('Angle of Attack α [deg]')
    plt.ylabel('Reynolds Number')
    plt.title(f'{quantity} Contour Plot')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{quantity}_contour.png')
    plt.show()

#plot_contour(dae51, quantity='CL')
#plot_contour(GEMINI, quantity='CL')
#plot_contour(ge51_25, quantity='CL')
#plot_contour(ge51_50, quantity='CL')
#plot_contour(ge51_75, quantity='CL')

