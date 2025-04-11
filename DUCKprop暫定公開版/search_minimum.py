# coding:utf-8

#改修3.ポーラーカーブの解析を実装。
#このコードを8割方実装し終えるまで、既存の翼弦自由のプログラム(interp2dのため動かず)が存在することを知らなかったためほぼ新造。

import numpy as np
import matplotlib.pyplot as plt
from load_and_interpolate import *
from scipy.interpolate import interp1d
from scipy.optimize import minimize

# Cl×Re＝Aの初期値、最大値、交差を設定
A_min = 0.04 #0.02以上の方が良い
A_max = 0.88 #少なくともClの最大値と規格化したReの最大値の乗算以下の値に
A_step = 0.07 #適切な値に

# ↑Aの項数が多いと、探索に非常に時間がかかる。項数は10程度が望ましい。項数が20近くになると反復(Iter)が100回を超える。

# グラフの背景（等高線図）の作成
    # メッシュグリッドの作成
Re_vals, CL_vals = np.meshgrid(dae51.Re, dae51.CL)

    # Cd × Re の等高線
Cd_Re_vals = dae51.CD_ReCL(Re_vals, CL_vals) * Re_vals

    # Cl × Re の等高線
Cl_Re_vals = CL_vals * Re_vals

    # 等高線図の作成
plt.figure(figsize=(8, 6))

    # Cd × Re の等高線の描画
cd_contour = plt.contourf(Re_vals, CL_vals, Cd_Re_vals, levels=50, cmap='jet', alpha=0.5)
plt.colorbar(cd_contour, label='CD × Re')

print()
print(f"start searching optimal points")

W = 100  # 重み

A_vals = np.arange(A_min, A_max, A_step)
n = len(A_vals)

# 初期値ベクトル（各点の初期[Re_i, Cl_i]）（テキトーに決めた）
x0 = []
for A in A_vals:
    plt.contour(Re_vals, CL_vals, Cl_Re_vals, levels=[A], colors='red', linewidths=0.5, linestyles='dashed')
    Re0 = A ** 0.5
    Cl0 = A / Re0 
    x0.extend([Re0, Cl0])
x0 = np.array(x0)

# 目的関数
def global_objective(x):
    total_cd_re = 0.0
    smoothness_penalty = 0.0

    for i in range(n):
        Re = x[2*i]
        Cl = x[2*i + 1]
        total_cd_re += dae51.CD_ReCL(Re, Cl) * Re

    # 滑らかさペナルティ (sin²θ)
    for i in range(1, n-1):
        x_prev = np.array([x[2*(i-1)], x[2*(i-1)+1]])
        x_curr = np.array([x[2*i], x[2*i+1]])
        x_next = np.array([x[2*(i+1)], x[2*(i+1)+1]])

        v1 = x_curr - x_prev
        v2 = x_next - x_curr

        norm_product = np.linalg.norm(v1) * np.linalg.norm(v2)
        if norm_product > 0:
            cos_theta = np.dot(v1, v2) / norm_product
            sin2_theta = 1 - cos_theta**2
            smoothness_penalty += sin2_theta

    # 第一項のReを非規格化するとバグる(SLSQP内部の問題？)ので第二項を規格化する
    return total_cd_re + W * smoothness_penalty / nc_Re 

# 制約条件：Re_i * Cl_i = A_i
constraints = []
for i in range(n):
    def make_constraint(i):
        return lambda x: x[2*i] * x[2*i + 1] - A_vals[i]
    constraints.append({'type': 'eq', 'fun': make_constraint(i)})

# 変数の範囲指定
bounds = []
for _ in range(n):
    bounds.extend([(0, 1), (0, 1.2)])  # Re, Cl

# 目的関数の履歴を記録するリスト
obj_history = []

# 進捗を表示するための関数
def print_progress(xk):
    obj_val = global_objective(xk)
    obj_history.append(obj_val)
    print(f"Iter {len(obj_history)}: Objective = {obj_val:.3f}")

# 最適化実行
result = minimize(
    global_objective, 
    x0, 
    method='SLSQP', 
    constraints=constraints, 
    bounds=bounds, 
    options={'disp': True, 'maxiter': 1000},
    callback=print_progress #進捗を表示
    )


# 最適化結果の処理
optimal_points = np.array(result.x).reshape(-1, 2)
optimal_cd_values = [dae51.CD_ReCL(Re, Cl) * Re * nc_Re for Re, Cl in optimal_points] #最適点のCd*Re（非規格化）


# 最適点の軌跡をプロット
if len(optimal_points) > 0:
    plt.plot(optimal_points[:, 0], optimal_points[:, 1], 'r-', linewidth=1, label="Optimal CD Path")

# グラフ設定
plt.xlabel('Reynolds Number (Re)')
plt.ylabel('Lift Coefficient (CL)')
plt.title('Contour Plot of CD × Re and Optimal Min CD Path')
plt.legend()
plt.grid(True)
plt.savefig('Contour Plot of CD × Re and Optimal Min CD Path.png')
plt.show()




#あるA（＝Cl×Re）におけるClを返す関数
# 各点の Cl × Re を計算
cl_re_values = optimal_points[:, 0] * optimal_points[:, 1]  # Re × Cl
cl_values = optimal_points[:, 1]  # Cl

# データを Cl × Re の昇順にソート（補間のため）
sorted_indices = np.argsort(cl_re_values)
cl_re_values = cl_re_values[sorted_indices]
cl_values = cl_values[sorted_indices]

# 補間の関数を作成（スプライン補間だとうまくいかなかった）
cl_int = interp1d(cl_re_values, cl_values, kind='linear', fill_value="extrapolate")

# スプライン関数を用いて補間曲線をプロット
cl_re_fine = np.linspace(min(cl_re_values), max(cl_re_values), 200)  # 滑らかな範囲
cl_fine = cl_int(cl_re_fine)

# グラフの描画
plt.figure(figsize=(8, 6))
plt.plot(cl_re_values, cl_values, 'ro', markersize=1, label="Optimal Cl Data")  # 元データ
plt.plot(cl_re_fine, cl_fine, 'b-', label="Linear Interpolation")  # 線形補間曲線

# 軸ラベルとタイトル
plt.xlabel('Cl × Re')
plt.ylabel('Cl')
plt.title('Spline Interpolation of Cl as a Function of Cl × Re')
plt.legend()
plt.grid(True)

# 画像として保存(↓のCd*Reのグラフと重ねてもいいですが、Cd*Reは非規格化していて桁があまりにも違うので分けています)
plt.savefig('optimal_cl_vs_cl_re.png')
plt.show()

#まとめ
optimaized_cl = lambda A: cl_spline(A)



# あるA（＝Cl×Re）におけるCd×Reを返す関数
optimal_cd_re_func = interp1d(A_vals, optimal_cd_values, kind='linear', fill_value="extrapolate")

# 補間曲線を描画
A_fine = np.linspace(min(A_vals), max(A_vals), 500)
optimal_cd_re_fine = optimal_cd_re_func(A_fine)
plt.figure(figsize=(8, 6))
plt.plot(A_vals, optimal_cd_values, 'bo', markersize=1, label="Optimal Cd × Re (Original)")
plt.plot(A_fine, optimal_cd_re_fine, 'r-', linewidth=1, label="Interpolated Function")

# グラフの設定
plt.xlabel("A = Cl × Re")
plt.ylabel("Optimal Cd × Re")
plt.title("Optimal Cd × Re vs A (Interpolated)")
plt.legend()
plt.grid(True)

# 画像として保存
plt.savefig("optimal_cd_re_vs_cl_re.png")    
plt.show()

#まとめ
optimaized_CDRe = lambda A: optimal_cd_re_func(A)