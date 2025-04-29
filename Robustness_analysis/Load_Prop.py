# coding:utf-8

import os
import re

import numpy
from scipy import interpolate
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

numpy.set_printoptions(
    threshold = numpy.inf,
    linewidth = 1e3
)

#プロペラ分割数
N = 50

# プロペラデータを置くディレクトリ
PROP_DIR = './prop'

def load_prop(file_name):
    print()
    # ファイルの検索
    file   = os.listdir('./'+PROP_DIR)
    for filename in file:
        if file_name not in filename:
                print('the data of different data included')
                print('check the contents of the folder and try again')
                exit()
        file = ['./'+PROP_DIR+'/'+filename for filename in file]
        print(f'start loading {file_name}...')

        large_table = None

    # データの読み込み
    for row_n, filename in enumerate(file):
        table = numpy.loadtxt(filename, skiprows = 1, usecols = [0, 1, 2], delimiter=',')
        # r, chord, theta の形の配列にする

        #読み込んだプロペラの図示
        fig=plt.figure(figsize=(12,5))

        plt.subplot(1,2,1)      
        plt.plot(table[:,0],table[:,1])
        plt.ylabel('Chord [mm]')
        plt.xlabel('r(distance from center of rotation)[m]')

        plt.subplot(1,2,2)
        plt.plot(table[:,0],table[:,2])
        plt.ylabel('theta [deg]')
        plt.xlabel('r(distance from center of rotation)[m]')

        #解析用データの準備：翼素に分割
        delta_B = ( table[-1,0] - table[0,0] ) / N

        #翼素の中心の回転半径
        r = []
        for i in range( N ):
            r.append( table[0,0] + delta_B * (i + 1/2))

        #翼素の幅を等しくするためにcsvのデジタルデータを補間してアナログ化
        interp_chord = interp1d( table[:,0], table[:,1], kind='linear', fill_value="extrapolate")
        interp_theta = interp1d( table[:,0], table[:,2], kind='linear', fill_value="extrapolate")
            
        #各rにおける翼弦長、取り付け角
        chord = interp_chord( r )
        theta = interp_theta( r )

        #分割したプロペラを重ねて図示（二つの曲線がそれぞれ重なっているか確認すること。両端は気にしなくてよろしい。）
        plt.subplot(1,2,1)      
        plt.plot( r , chord )

        plt.subplot(1,2,2)
        plt.plot( r , theta )

        fig.savefig('Loaded_prop.png')
        plt.show()

        large_table = numpy.column_stack([ r , chord , theta ])

        #各翼素の中心の位置ベクトル
        CP = numpy.array([
            numpy.array([0, r[i], 0]) for i in range(N)
        ])

    print('complete')
    print()
    return large_table , delta_B , CP

large_table , delta_B , CP = load_prop('for_design_at2025-04-08_160014') #ファイル名をそのまま貼り付け