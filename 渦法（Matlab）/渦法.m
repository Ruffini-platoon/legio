%% constant_and_setting

% 物理定数(@ Altitude 0m, Temp 30度)
RHO = 1.164;
NU  = 1.60e-5;

% 設計依存の定数
B     = 2; % ブレード数
U_inf = 7.9; % 機速[m/s]
RPM   = 143; % ペラ回転数[rpm]
%Power = ; % 出力[W]
Thrust = 21; % 推力[N]
Omega = deg2rad(RPM*360) / 60; % 回転数[rad/sec]
R_out  = 1.48; % ペラ端回転半径[m]
R_in   = 0.059; % ペラ根回転半径[m]


% 解析のための定数
N  = 100; % ペラ分割数
M  = 100; % 馬蹄渦分割
%dt = 0.05; % 馬蹄渦分割時間幅(廃止)

%非一様な馬蹄渦分割時間幅
t1 = 1/2000;                             %ペラ付近のt_step
t2 = 1/20;                               %離れたところのt_step
lambda = 40;                             %収束の速さのパラメータ
t = [0:M];
t_step = 1./(1+exp(t/lambda));
t_step = 1-2*t_step;
t_step = (t_step)*(t2-t1)+t1;
t = cumsum(t_step);
for i = 1:M+1
    t(i) = t(i)-t_step(i);
end

function [CP,delta_B,r] = sep_blade(R_out, R_in, N)
    delta_B = (R_out - R_in) / N;
    r = linspace(R_in ,R_out-delta_B ,N) + delta_B/2;
    r = r(:);
    CP = [zeros(N,1), r, zeros(N,1)];
end

function DP = sep_vortex(delta_B, r, U_inf, Omega, B, N, M, t)
    r = r - delta_B/2;
    r = [r; r(end) + delta_B];

    DP = zeros(B, N+1, M+1, 3);
    for bl = 1:B
        for j = 1:(N+1)
            for k = 1:(M+1)
                %t = (k-1)*dt;
                theta = Omega * t(k) + 2*pi*(bl-1)/B;

                x = -U_inf * t(k);
                y = r(j) * cos(-theta);
                z = r(j) * sin(-theta);
                DP(bl, j, k, :) = [x, y, z];
            end
        end
    end
end

%代表点の座標
[CP, delta_B, r]= sep_blade(R_out, R_in, N);
DP = sep_vortex(delta_B, r, U_inf, Omega, B, N, M, t);

%プロペラ後流の長さ
Length_of_PropWake = U_inf*t(end) %ペラ半径の10倍以上の値になっているか確認

%% load_and_interpolate

function airfoil = AirFoil(airfoil_name)
    % AIRFOIL_DIR はグローバル定数
    AIRFOIL_DIR = './airfoil';

    % ===== ファイルの検索 =====
    folder_path = fullfile(AIRFOIL_DIR, airfoil_name);
    files = dir(fullfile(folder_path, '*'));
    files = files(~[files.isdir]); % フォルダを除外

    % ファイル名チェック
    for i = 1:length(files)
        if ~contains(files(i).name, [airfoil_name '_T1'])
            error('the data of different airfoil included\ncheck the contents of the folder and try again');
        end
    end

    fprintf('start loading %s...\n', airfoil_name);

    large_table2 = [];

    figure;
    % ===== データ読み込み =====
    for i = 1:length(files)
        filename = fullfile(folder_path, files(i).name);

        % ファイル名から Re 数を抽出
        token = regexp(files(i).name, 'Re\d+\.\d+', 'match');
        if isempty(token)
            error('Re number not found in file name: %s', files(i).name);
        end
        Re_val = str2double(token{1}(3:end)) * 1e6;
        fprintf('    loading file of Re=%g\n', Re_val);

        % データ読み込み (alpha, CL, CD)
        table_data = readmatrix(filename, 'NumHeaderLines', 10);
        table_data = table_data(:, 1:3); % α, CL, CD

        % Re 列を追加
        Re_col = repmat(Re_val, size(table_data, 1), 1);
        table_data = [Re_col, table_data]; % Re, α, CL, CD

        % グラフ描画（CD vs CL）
        plot(table_data(:, 4), table_data(:, 3)); hold on;

        % large_table2 に縦結合
        large_table2 = [large_table2; table_data];
    end

    saveas(gcf, 'polar.png');
    fprintf('complete\n');

    % グリッド作成
    nc_Re = 3e5; %規格化定数

    Re_grid = [1e4 / nc_Re, 3e5 / nc_Re, 500];
    CL_grid = [0, 1.2, 500];

    Re_vec = linspace(Re_grid(1), Re_grid(2), Re_grid(3));
    CL_vec = linspace(CL_grid(1), CL_grid(2), CL_grid(3));

    [Re_mesh, CL_mesh] = meshgrid(Re_vec, CL_vec);

    % alphaの補間データ作成
    alpha_grid = griddata( ...
        large_table2(:, 1) / nc_Re, large_table2(:, 3), large_table2(:, 2), ...
        Re_mesh, CL_mesh, 'cubic' ...
    );

    % CDの補間データ作成
    CD_grid = griddata( ...
        large_table2(:, 1) / nc_Re, large_table2(:, 3), large_table2(:, 4), ...
        Re_mesh, CL_mesh, 'cubic' ...
    );

    % griddedInterpolantで線形補間＋線形外挿
    F_alpha = griddedInterpolant({Re_vec, CL_vec}, alpha_grid', 'cubic', 'linear');
    F_CD = griddedInterpolant({Re_vec, CL_vec}, CD_grid', 'cubic', 'linear');

    airfoil.alpha_ReCL = @(Re, CL) F_alpha(Re / nc_Re, CL);%規格化していないReを引数に
    airfoil.CD_ReCL = @(Re, CL) F_CD(Re, CL);%規格化したReを引数に

    % 保存
    airfoil.Re = Re_vec;
    airfoil.CL = CL_vec;

    fprintf('airfoil data successfully loaded\n');
end

dae51 = AirFoil('DAE51');

%Re数を規格化しているので注意‼

%dae51.CD_ReCL(0.5,0.05) %確認用

%% search_optimal_path

%Cl×Re＝Aの初期値、最大値、交差を設定
A_min = 0.01; 
A_max = 0.86; %少なくともClの最大値と規格化したReの最大値の乗算以下の値に
A_step = 0.05; %適切な値に。Aの項数が多いと、探索に非常に時間がかかる。

% グラフの背景（等高線図）の作成
    % メッシュグリッドの作成
[Re_vals, CL_vals] = meshgrid(dae51.Re, dae51.CL);

    % Cd × Re の等高線
Cd_Re_vals = dae51.CD_ReCL(Re_vals, CL_vals) .* Re_vals;

    % Cl × Re の等高線
Cl_Re_vals = CL_vals .* Re_vals;

figure('Position',[100 100 800 600]); % 図のサイズ指定

% --- 背景 Cd×Re 等高線 ---
ax1 = axes;
contourf(ax1, Re_vals, CL_vals, Cd_Re_vals, 100, 'LineColor','none'); 
colormap(ax1, jet);
colorbar; ylabel(colorbar,'CD × Re');  
grid on; hold on;

disp(' ');
disp('start searching optimal points');

% Cl×Re = A の等高線
A_vals = A_min : A_step : A_max;
%A_vals = [0.02,0.03,0.04,0.05,0.06,0.07,0.09,0.12,0.16,0.19,0.37,0.45,0.51,0.57,0.66,0.75,0.84];
n = length(A_vals);
W = 300; % 重み
x0 = []; % 初期値ベクトル

% --- Cl×Re 等高線用 axes（透明背景） ---
ax2 = axes('Position', ax1.Position);
hold(ax2,'on');

nc_Re = 3e5; %関数内で登場させてしまったので再掲

for i = 1:length(A_vals)
    A = A_vals(i);
    
    % Cl*Re = A の等高線を描く
    contour(ax2, Re_vals, CL_vals, Cl_Re_vals, [A A], 'LineColor', 'r', ...
        'LineWidth', 0.5, 'LineStyle', '--');
    
    % 初期値 Re0, Cl0
    Re0 = sqrt(A);
    Cl0 = A / Re0;
    
    % x0 に追加（行ベクトルとして）
    x0 = [x0, Re0, Cl0];
end

% --- ax2 の設定 ---
ax2.Color = 'none'; % 背景透明
ax2.XLim = ax1.XLim; ax2.YLim = ax1.YLim; % 軸範囲同期
ax2.XLabel.String = 'Re'; ax2.YLabel.String = 'CL';
ax2.Box = 'on';

%目的関数
function f = global_objective(x, dae51, n, W, nc_Re)
    % x: [Re1, Cl1, Re2, Cl2, ..., Re_n, Cl_n]
    % dae51: AirFoil 構造体
    % n: 等高線の数
    % W: 滑らかさペナルティの重み
    % nc_Re: Re の規格化定数

    total_cd_re = 0;
    smoothness_penalty = 0;

    % --- total_cd_re の計算 ---
    for i = 1:n
        Re = x(2*i - 1);
        Cl = x(2*i);
        total_cd_re = total_cd_re + dae51.CD_ReCL(Re, Cl) * Re;
    end

    % --- 滑らかさペナルティ (sin^2 θ) ---
    for i = 2:(n-1)
        x_prev = [x(2*(i-1)-1), x(2*(i-1))];
        x_curr = [x(2*i - 1), x(2*i)];
        x_next = [x(2*(i+1)-1), x(2*(i+1))];

        v1 = x_curr - x_prev;
        v2 = x_next - x_curr;

        norm_product = norm(v1) * norm(v2);
        if norm_product > 0
            cos_theta = dot(v1, v2) / norm_product;
            sin2_theta = 1 - cos_theta^2;
            smoothness_penalty = smoothness_penalty + sin2_theta;
        end
    end
    f = total_cd_re + W * smoothness_penalty/nc_Re;
end

%制約条件
function [c, ceq] = constraint_eq(x, A_vals, n)
    % c: 不等式制約 (c <= 0) -> ここでは無し
    % ceq: 等式制約 (ceq = 0)

    c = []; % 不等式は無し
    ceq = zeros(1, n); % n 個の等式制約

    for i = 1:n
        Re = x(2*i - 1);
        Cl = x(2*i);
        ceq(i) = Re * Cl - A_vals(i);
    end
end

% 変数の範囲指定（下限と上限）
lb = repmat([1/30, 0], 1, n);    % Re >= 1/30, Cl >= 0
ub = repmat([1, 1.1], 1, n);     % Re <= 1, Cl <= 1.1

% 最適化オプション
options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...           % SLSQP 相当
    'Display','iter', ...
    'MaxIterations', 1000, ...       % 反復制限
    'MaxFunctionEvaluations', 1e10); % 関数評価制限

% 最適化実行
[x_opt, fval] = fmincon(@(x) global_objective(x, dae51, n, W, nc_Re), ...
                        x0, ...   % 初期値
                        [], [], [], [], ... % 線形制約なし
                        lb, ub, ...         % 変数の範囲
                        @(x) constraint_eq(x, A_vals, n), ... % 非線形制約
                        options);

% 最適化結果の処理
optimal_points = reshape(x_opt, 2, []).';  % n×2 の行列に変換: [Re_i, Cl_i]

% 最適点の Cd*Re（非規格化）
optimal_cd_values = zeros(size(optimal_points,1),1);
for i = 1:size(optimal_points,1)
    Re = optimal_points(i,1);
    Cl = optimal_points(i,2);
    optimal_cd_values(i) = dae51.CD_ReCL(Re, Cl) * Re * nc_Re; % 非規格化
end

% 背景の等高線は前に描いた Cd*Re contourf がある想定
if ~isempty(optimal_points)
    plot(optimal_points(:,1), optimal_points(:,2), 'r-', 'LineWidth', 1, 'DisplayName', 'Optimal CD Path');
end

% グラフ設定
xlabel('Reynolds Number (Re)');
ylabel('Lift Coefficient (CL)');
title(sprintf('Optimal Min CD Path  W = %.1f  A:(%.2f , %.2f , %.2f)', W, A_min, A_max, A_step));
grid on;

% 保存
saveas(gcf, 'Contour Plot of CD × Re and Optimal Min CD Path.png');

%結果表示
% disp('Optimal points [Re, CL]:');
% disp(optimal_points);
% 
% disp('Optimal Cd*Re values:');
% disp(optimal_cd_values);

% --- Cl × Re = A における Cl 補間関数作成 ---

% 各点の Cl × Re を計算
cl_re_values = optimal_points(:,1) .* optimal_points(:,2);  % Re × Cl
cl_values    = optimal_points(:,2);                         % Cl

% Cl × Re の昇順にソート（補間用）
[cl_re_values, sorted_indices] = sort(cl_re_values);
cl_values = cl_values(sorted_indices);

% 線形補間関数を作成
cl_int = @(A) interp1(cl_re_values, cl_values, A, 'linear', 'extrap');

% 補間曲線を描画
cl_re_fine = linspace(min(cl_re_values), max(cl_re_values), 200);
cl_fine = cl_int(cl_re_fine);

figure;
plot(cl_re_values, cl_values, 'ro', 'MarkerSize', 3, 'DisplayName', 'Optimal Cl Data'); hold on;
plot(cl_re_fine, cl_fine, 'b-', 'LineWidth', 1, 'DisplayName', 'Linear Interpolation');
xlabel('Cl × Re'); ylabel('Cl');
title('Linear Interpolation of Cl as a Function of Cl × Re');
legend('show'); grid on;

% --- 関数としてまとめる ---
optimaized_cl = cl_int;  % @(A) interp1(...)

% --- Cl × Re = A における Cd × Re 補間関数作成 ---

% 線形補間関数を作成
optimal_cd_re_func = @(A) interp1(A_vals, optimal_cd_values, A, 'linear', 'extrap');

% 補間曲線を描画
A_fine = linspace(min(A_vals), max(A_vals), 500);
optimal_cd_re_fine = optimal_cd_re_func(A_fine);

figure;
plot(A_vals, optimal_cd_values, 'bo', 'MarkerSize', 3, 'DisplayName', 'Optimal Cd × Re (Original)'); hold on;
plot(A_fine, optimal_cd_re_fine, 'r-', 'LineWidth', 1, 'DisplayName', 'Interpolated Function');
xlabel('A = Cl × Re'); ylabel('Optimal Cd × Re');
title('Optimal Cd × Re vs A (Interpolated)');
legend('show'); grid on;

% --- 関数としてまとめる ---
optimaized_CDRe = optimal_cd_re_func;

%% represent_with_gamma

% 馬蹄渦が引き起こす速度計算関数 f
function fval = f(a,b)
        l = b - a;
        al = cross(a,l,2);
        norm_al = vecnorm(al, 2, 2);
        e_a = a ./ vecnorm(a, 2, 2);
        e_b = b ./ vecnorm(b, 2, 2);
        fval = zeros(size(a));
        zero_mask = (norm_al == 0);
        denom = norm_al.^2;
        numerator = al ./ denom;
        dot_val = dot((e_b - e_a), l, 2);
        fval(~zero_mask,:) = (1/(4*pi)) * numerator(~zero_mask,:) .* dot_val(~zero_mask,:);
end

% 誘導速度行列 X,Y,Z の計算
disp('')
fprintf('start computing biot-savart coefficient matrix\n');

% 配列生成
a = zeros(B, N, N, M, 3);
b = zeros(B, N, N, M, 3);
a_dash = zeros(B, N, N, M, 3);
b_dash = zeros(B, N, N, M, 3);

for bl = 1:B
    for i = 1:N
        for j = 1:N
            for k = 1:M
                a(bl,i,j,k,:) = squeeze(DP(bl,j,k,:))' - CP(i,:);
                b(bl,i,j,k,:) = squeeze(DP(bl,j,k+1,:))' - CP(i,:);
                a_dash(bl,i,j,k,:) = squeeze(DP(bl,j+1,k,:))' - CP(i,:);
                b_dash(bl,i,j,k,:) = squeeze(DP(bl,j+1,k+1,:))' - CP(i,:);
            end
        end
    end
end
fprintf('  complete setting relative position\n');

X = zeros(N,N);
Y = zeros(N,N);
Z = zeros(N,N);

fprintf('  computing coefficient matrix\n');
for i = 1:N
    fprintf('    %.2f%%...\n', 100*i/N);
    for j = 1:N
        s1 = zeros(1,3);
        for bl = 1:B
            for k = 1:M
                s1 = s1 + f(squeeze(a(bl,i,j,k,:))', squeeze(b(bl,i,j,k,:))');
            end
        end
        s2 = zeros(1,3);
        for bl = 1:B
            s2 = s2 + f(squeeze(a(bl,i,j,1,:))', squeeze(a_dash(bl,i,j,1,:))');
        end

        s3 = zeros(1,3);
        for bl = 1:B
            for k = 1:M
                s3 = s3 + f(squeeze(a_dash(bl,i,j,k,:))', squeeze(b_dash(bl,i,j,k,:))');
            end
        end

        val = -s1 + s2 + s3;
        X(i,j) = val(1);
        Y(i,j) = val(2);
        Z(i,j) = val(3);
    end
end
fprintf('  complete\n');

u_ind = @(gamma) X * gamma(:);

w_ind = @(gamma) Z * gamma(:);

%誘導速度の大きさ
F_IV = @(gamma) sqrt(u_ind(gamma).^2 + w_ind(gamma).^2);

% 対気速度の x 成分 (式6)
F_Up = @(gamma) U_inf - u_ind(gamma);

% 対気速度の z 成分 (式5)
F_Ut = @(gamma) Omega * r - w_ind(gamma);

% 対気速度の大きさ (式7, 18)
F_V  = @(gamma) sqrt(F_Up(gamma).^2 + F_Ut(gamma).^2);

% 誘導迎え角 (式8, 19)
F_phi = @(gamma) atan(F_Up(gamma) ./ F_Ut(gamma));


% 局所揚力 (式9, 20)
F_L = @(gamma) delta_B * RHO * (F_V(gamma) .* gamma(:));

% Re*Cl (式39)
F_ReCl = @(gamma) (2 .* F_L(gamma)) ./ (RHO .* NU .* F_V(gamma) .* delta_B .* nc_Re);

% 揚力係数 (式14, 21, 44)
F_CL = @(gamma) optimaized_cl(F_ReCl(gamma));

% レイノルズ数 (式45)
F_Re = @(gamma) (F_ReCl(gamma) * nc_Re) ./ F_CL(gamma);

% 抗力係数 (式15, 22, 47)
F_CD = @(gamma) optimaized_CDRe(F_ReCl(gamma)) ./ F_Re(gamma);

% 翼弦長 (式16, 23, 48)
F_C = @(gamma) (NU .* F_ReCl(gamma) .* nc_Re) ./ (F_V(gamma) .* F_CL(gamma));

% 迎角 (式17, 24, 49)
F_alpha = @(gamma) dae51.alpha_ReCL(F_Re(gamma), F_CL(gamma));

% gamma のみの関数としての dD (式10, 25)
F_D = @(gamma) 0.5 * RHO * (F_V(gamma).^2) .* F_C(gamma) .* delta_B .* F_CD(gamma);

% 推力 (式12, 27)
F_T = @(gamma) B * sum(F_L(gamma) .* cos(F_phi(gamma)) - F_D(gamma) .* sin(F_phi(gamma)));

minusT = @(gamma) -F_T(gamma);

% パワー (式13, 28)
F_P = @(gamma) B * Omega * sum((F_D(gamma) .* cos(F_phi(gamma)) + F_L(gamma) .* sin(F_phi(gamma))) .* r);

%% main

gamma0 = sin(linspace(0,1,N) * pi) + 0.01; %循環は行ベクトル

fprintf('start optimizing gamma...\n');

% 下限・上限 (いらないかも)
lb = zeros(N,1);
ub = 2 * ones(N,1);

% 不等式制約
function c = ineqcon_gamma(gamma)

    % (1) 二階差分: 循環分布が凸になるように
    k = 1:length(gamma)-2; 
    con1 = -2*gamma(k+1) + gamma(k) + gamma(k+2);
    % (2) 翼根循環の下限制約
    %constraint = 0.01;
    %con2 = constraint - gamma(1);

    % まとめ (c <= 0)
    %c = [con1(:); con2];
    c = [con1(:)];
end

% 等式制約
eqcon = @(gamma) F_T(gamma) - Thrust;

nonlcon = @(gamma) deal(ineqcon_gamma(gamma), eqcon(gamma));

% オプション
options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','iter', ...
    'MaxIterations', 1e8, ...
    'OptimalityTolerance', 1e-8, ...
    'MaxFunctionEvaluations', 1e10);

% 最適化実行
[x_opt, fval, exitflag, output] = fmincon(@(x) F_P(x), gamma0, ...
    [], [], [], [], lb, ub, nonlcon, options);

fprintf('optimization completed\n');
gamma = x_opt;  % 最適解

% --- 結果計算 ---
CL = F_CL(gamma);
CD = F_CD(gamma);
V  = F_V(gamma);
dL = F_L(gamma);
dD = F_D(gamma);
Re = (F_V(gamma).*F_C(gamma))/NU;
phi = rad2deg(F_phi(gamma)); % 誘導迎え角 
alpha = F_alpha(gamma);      % 有効迎え角 
theta = phi + alpha;         % 取付迎角
thrust = F_T(gamma);
power = F_P(gamma);
eta = thrust * U_inf / power;
chord = F_C(gamma);

fprintf('thrust: %g\n', thrust);
fprintf('power : %g\n', power);
fprintf('eta   : %g\n', eta);

% --- テーブル作成（行列: 行 = 各翼素） ---
% 列順: r, chord, theta, phi, alpha, gamma, Re, CL, CD, V, dL, dD
T = [r(:), chord(:), theta(:), phi(:), alpha(:), gamma(:), Re(:), CL(:), CD(:), V(:), dL(:), dD(:)];

% ヘッダ文字列
nowstr = datestr(now,'yyyy-mm-dd_HHMMSS');
folder = ['prop_result_' nowstr];
if ~exist(folder,'dir'), mkdir(folder); end

txtfile = fullfile(folder, ['result_at_' nowstr '.txt']);
%{
hdr = sprintf(['RHO   : %g\nNU    : %g\nB     : %d\nU_inf : %g\nRPM   : %g\npower : %g\nR_out : %g\nR_in  : %g\n' ...
    'N     : %d\nM     : %d\ndt    : %g\nthrust: %g\neta   : %g\noptsuccess:%s\n'], ...
    RHO, NU, B, U_inf, RPM, power, R_out, R_in, N, M, dt, thrust, eta, output.message);
%}
hdr = sprintf(['RHO   : %g\nNU    : %g\nB     : %d\nU_inf : %g\nRPM   : %g\npower : %g\nR_out : %g\nR_in  : %g\n' ...
    'N     : %d\nM     : %d\nt1    : %d\nt2    : %g\nthrust: %g\neta   : %g\noptsuccess:%s\n'], ...
    RHO, NU, B, U_inf, RPM, power, R_out, R_in, N, M, t1, t2, thrust, eta, output.message);
colnames = '半径r[m],翼弦長chord[m],取付角theta[deg],流入角phi[deg],迎角alpha[deg],循環gamma,Re,CL,CD,V,dL,dD';

% 保存（数値はCSVで、ヘッダは先頭に）
fid = fopen(txtfile,'w');
fprintf(fid, '%s\n', hdr);
fprintf(fid, '%s\n', colnames);
fclose(fid);
writematrix(T, txtfile, 'Delimiter',',','WriteMode','append');

% --- 設計用CSV (補間) ---
r_spline = (0:0.001:R_out).';
chord_spline_list = interp1(r, chord, r_spline, 'linear', 'extrap');
theta_spline_list = interp1(r, theta, r_spline, 'linear', 'extrap');

csvfile = fullfile(folder, ['for_design_at' nowstr '.csv']);
writematrix([r_spline, chord_spline_list, theta_spline_list], csvfile);

% --- 図1: 3x3 表示 ---
fig1 = figure('Visible','off','Position',[100 100 1200 840]);
subplot(3,3,1); plot(r, chord); title('Chord(chord length)');
subplot(3,3,2); plot(r, theta); title('Theta(mounting angle)');
subplot(3,3,3); plot(r, gamma); title('gamma(circulation)');
subplot(3,3,4); plot(r, Re); title('Re');
subplot(3,3,5); plot(r, CL); title('CL');
subplot(3,3,6); plot(r, CD); title('CD');
subplot(3,3,7); plot(r, V); title('V');
subplot(3,3,8); plot(r, dL); title('dL');
subplot(3,3,9); plot(r, dD); title('dD');
saveas(fig1, fullfile(folder, [nowstr '_propfig.png']));
close(fig1);

% --- 図2: thesis style 4x2 ---
fig2 = figure('Visible','off','Position',[100 100 1000 1300]);
tiledlayout(4,2);
nexttile; plot(r,chord); ylabel('Chord [m]');
nexttile; plot(r,theta); ylabel('Theta [deg]');
nexttile; plot(r,gamma); ylabel('gamma [m^2/s]');
nexttile; plot(r,Re); ylabel('Re');
nexttile; plot(r,CL); ylabel('CL');
nexttile; plot(r,CD); ylabel('CD');
nexttile; plot(r,V); ylabel('V [m/s]'); xlabel('r [m]');
nexttile; plot(r,alpha); ylabel('alpha [deg]'); xlabel('r [m]');
saveas(fig2, fullfile(folder, [nowstr '_propfig_thesis.png']));
close(fig2);