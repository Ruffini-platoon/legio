%% Load_Prop.pyの内容

% 解析したいプロペラのデータを同じフォルダに入れ、以下にそのファイル名を入力
filename = '26代ペラ.csv';
Prop = readtable(filename);
% disp(Prop);%確認用

% 分割数(任意)
N = 100;

% tableの半径の値が翼素の中心の回転半径ならば以下の式
pre_delta_B = (Prop{height(Prop),1}-Prop{1,1})/(height(Prop)-1);
delta_B = (Prop{height(Prop),1}-Prop{1,1}+pre_delta_B)/N;
% リブの回転半径（翼素の端の回転半径）ならば以下の式を使う。
% delta_B = (Prop{200,1}-Prop{1,1})/N;

% 翼素の中心の回転半径（tableの半径の値がリブの回転半径ならば適宜変更すること）
r = linspace(Prop{1,1}-0.5*(pre_delta_B-delta_B),Prop{height(Prop),1}+0.5*(pre_delta_B-delta_B),N);

interp_chord = griddedInterpolant(Prop{:,1}, Prop{:,2}, 'linear', 'linear');
interp_theta = griddedInterpolant(Prop{:,1}, Prop{:,3}, 'linear','linear');

chord = interp_chord(r);
theta = interp_theta(r);

r = r(:);         % 縦ベクトル化
chord = chord(:); % 縦ベクトル化
theta = theta(:); % 縦ベクトル化

large_table = [r, chord, theta]; % N×3行列

CP = [zeros(N,1), r, zeros(N,1)];% N×3の行列

%% load_and_interpolate.pyの内容。

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

        % ファイル名から Re 数を抽出（例: Re2.345 → 2.345×1e6）
        token = regexp(files(i).name, 'Re\d\.\d{3}', 'match');
        if isempty(token)
            error('Re number not found in file name: %s', files(i).name);
        end
        Re_val = str2double(token{1}(3:end)) * 1e6;
        fprintf('    loading file of Re=%g\n', Re_val);

        % データ読み込み (alpha, CL, CD)
        table_data = readmatrix(filename, 'NumHeaderLines', 11);%xflrから持ってきたファイルが新しいものなら10に、古いものなら11に
        table_data = table_data(:, 1:3); % α, CL, CD

        % Re 列を追加
        Re_col = repmat(Re_val, size(table_data, 1), 1);
        table_data = [Re_col, table_data]; % Re, α, CL, CD

        % グラフ描画（CD vs CL）
        plot(table_data(:, 4), table_data(:, 3)); hold on;

        % large_table に縦結合
        large_table2 = [large_table2; table_data];
    end

    saveas(gcf, 'polar.png');
    fprintf('complete\n');

    % グリッド作成
    Re_grid = [1e4, 3e5, 50];
    alpha_grid = [-10, 25, 50];

    Re_vec = linspace(Re_grid(1), Re_grid(2), Re_grid(3));
    alpha_vec = linspace(alpha_grid(1), alpha_grid(2), alpha_grid(3));

    [alpha_mesh, Re_mesh] = meshgrid(alpha_vec, Re_vec);

    % CLの補間データ作成
    CL_grid = griddata( ...
        large_table2(:, 2), large_table2(:, 1), large_table2(:, 3), ...
        alpha_mesh, Re_mesh, 'cubic' ...
    );

    % CDの補間データ作成
    CD_grid = griddata( ...
        large_table2(:, 2), large_table2(:, 1), large_table2(:, 4), ...
        alpha_mesh, Re_mesh, 'cubic' ...
    );

    % griddedInterpolantで線形補間＋線形外挿
    F_CL = griddedInterpolant({alpha_vec, Re_vec}, CL_grid', 'linear', 'linear');
    F_CD = griddedInterpolant({alpha_vec, Re_vec}, CD_grid', 'linear', 'linear');

    airfoil.CL_alphaRe = @(alpha, Re) F_CL(alpha, Re);
    airfoil.CD_alphaRe = @(alpha, Re) F_CD(alpha, Re);

    % 保存
    airfoil.Re = Re_vec;
    airfoil.alpha = alpha_vec;

    fprintf('airfoil data successfully loaded\n');
end

dae51 = AirFoil('DAE51');
%DAE51 = AirFoil('DAE51');
%ga43 = AirFoil('ga43');
%ga42 = AirFoil('ga42');
%ga41 = AirFoil('ga41');
%ga40 = AirFoil('ga40');
%ga39 = AirFoil('ga39');
%ga38 = AirFoil('ga38');
%ga37 = AirFoil('ga37');
%ga36 = AirFoil('ga36');
%ga35 = AirFoil('ga35');
%ga34 = AirFoil('ga34');
%ga33 = AirFoil('ga33');
%ga32 = AirFoil('ga32');
%ga31 = AirFoil('ga31');
%ga30 = AirFoil('ga30');
%ga29 = AirFoil('ga29');
%ga28 = AirFoil('ga28');
%ga27 = AirFoil('ga27');
%ga26 = AirFoil('ga26');
%ga25 = AirFoil('ga25');
%ga24 = AirFoil('ga24');
%ga23 = AirFoil('ga23');
%ga22 = AirFoil('ga22');
%ga21 = AirFoil('ga21');
%ga20 = AirFoil('ga20');
%ga19 = AirFoil('ga19');
%ga18 = AirFoil('ga18');
%ga17 = AirFoil('ga17');
%ga16 = AirFoil('ga16');
%ga15 = AirFoil('ga15');
%ga14 = AirFoil('ga14');
%ga13 = AirFoil('ga13');
%ga12 = AirFoil('ga12');
%ga11 = AirFoil('ga11');
%ga10 = AirFoil('ga10');
%ga9 = AirFoil('ga9');
%ga8 = AirFoil('ga8');
%ga7 = AirFoil('ga7');
%ga6 = AirFoil('ga6');
%ga5 = AirFoil('ga5');
%ga4 = AirFoil('ga4');
%ga3 = AirFoil('ga3');
%ga2 = AirFoil('ga2');
%ga1 = AirFoil('ga1');
%RootFOIL = AirFoil('RootFOIL');

%% Analysis.pyの内容

function [T, P, eta, norm_val] = analysis_at(U_inf, RPM, large_table, delta_B, CP, N, ...
    DAE51)
    % 物理定数
    RHO = 1.164;% @ Altitude 0m, Temp 30度
    NU = 1.60e-5;% @ Altitude 0m, Temp 30度
    B = 2;  % ブレード数
    M  = 100; % 馬蹄渦分割
    %dt = 0.05; % 馬蹄渦分割時間幅(廃止)
    eps = 0.05; %循環の滑らかさペナルティ閾値

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

    Omega = deg2rad(RPM*360) / 60;

    r = large_table(:,1);
    chord = large_table(:,2);
    theta = large_table(:,3);

    % 馬蹄渦の分割点生成
    function DP = sep_vortex(CP, delta_B, r, U_inf, Omega, B, N, M, t)
        r_shifted = r - delta_B/2;
        r_ext = [r_shifted; r_shifted(end) + delta_B];

        DP = zeros(B, N+1, M+1, 3);
        for bl = 1:B
            for j = 1:(N+1)
                for k = 1:(M+1)
                    %t = (k-1)*dt;
                    theta_ = 2*pi/B * (bl-1) + Omega*t(k);
                    x = - U_inf * t(k);
                    y = r_ext(j) * cos(-theta_);
                    z = r_ext(j) * sin(-theta_);
                    DP(bl,j,k,:) = [x, y, z];
                end
            end
        end
    end

    DP = sep_vortex(CP, delta_B, r, U_inf, Omega, B, N, M, t);

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
        fprintf('    %.2f%%...\n', floor(100*i/N));
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

    % u_ind, w_indの関数ハンドル定義
    u_ind = @(gamma) X * gamma(:);
    w_ind = @(gamma) Z * gamma(:);

    % 各関数定義
    F_Up = @(gamma) U_inf - u_ind(gamma);
    F_Ut = @(gamma) Omega * r - w_ind(gamma);
    F_V  = @(gamma) sqrt(F_Up(gamma).^2 + F_Ut(gamma).^2);
    F_phi = @(gamma) atan(F_Up(gamma) ./ F_Ut(gamma));
    F_alpha = @(gamma) theta - rad2deg(F_phi(gamma));
    F_Re = @(gamma) F_V(gamma) .* chord / NU;

    %単一翼型の場合
    F_CL = @(gamma) DAE51.CL_alphaRe(F_alpha(gamma), F_Re(gamma));
    F_CD = @(gamma) DAE51.CD_alphaRe(F_alpha(gamma), F_Re(gamma));

    %{
    %複数翼型の場合
    function result = F_CL(gamma,r)
        result = zeros(size(gamma));

        F_alpha_result = F_alpha(gamma);
        F_Re_result = F_Re(gamma);
        
        for k = 1:length(gamma)
            rad = r(k);     % rのk番目のスカラー値
        
        % 条件分岐
            if rad < 0.066
                result(k) = RootFOIL.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif (rad >= 0.066 && rad < 0.073)
                result(k) = ga43.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif (rad >= 0.073 && rad < 0.08)
                result(k) = ga42.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif (rad >= 0.08 && rad < 0.087)
                result(k) = ga41.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.087 && rad < 0.094
                result(k) = ga40.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.094 && rad < 0.101
                result(k) = ga39.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.101 && rad < 0.108
                result(k) = ga38.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.108 && rad < 0.115
                result(k) = ga37.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.115 && rad < 0.122
                result(k) = ga36.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.122 && rad < 0.129
                result(k) = ga35.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.129 && rad < 0.136
                result(k) = ga34.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.136 && rad < 0.143
                result(k) = ga33.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.143 && rad < 0.15
                result(k) = ga32.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.15 && rad < 0.157
                result(k) = ga31.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.157 && rad < 0.164
                result(k) = ga30.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.164 && rad < 0.172
                result(k) = ga29.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.172 && rad < 0.179
                result(k) = ga28.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.179 && rad < 0.186
                result(k) = ga27.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.186 && rad < 0.193
                result(k) = ga26.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.193 && rad < 0.2
                result(k) = ga25.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.2 && rad < 0.208
                result(k) = ga24.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.208 && rad < 0.215
                result(k) = ga23.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.215 && rad < 0.222
                result(k) = ga22.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.222 && rad < 0.229
                result(k) = ga21.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.229 && rad < 0.236
                result(k) = ga20.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.236 && rad < 0.244
                result(k) = ga19.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.244 && rad < 0.251
                result(k) = ga18.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.251 && rad < 0.258
                result(k) = ga17.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.258 && rad < 0.265
                result(k) = ga16.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.265 && rad < 0.272
                result(k) = ga15.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.272 && rad < 0.279
                result(k) = ga14.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.279 && rad < 0.286
                result(k) = ga13.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.286 && rad < 0.293
                result(k) = ga12.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.293 && rad < 0.3
                result(k) = ga11.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.3 && rad < 0.307
                result(k) = ga10.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.307 && rad < 0.314
                result(k) = ga9.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.314 && rad < 0.321
                result(k) = ga8.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.321 && rad < 0.329
                result(k) = ga7.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.329 && rad < 0.336
                result(k) = ga6.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.336 && rad < 0.343
                result(k) = ga5.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.343 && rad < 0.35
                result(k) = ga4.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.35 && rad < 0.357
                result(k) = ga3.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.357 && rad < 0.364
                result(k) = ga2.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.364 && rad < 0.371
                result(k) = ga1.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.371
                result(k) = DAE51.CL_alphaRe(F_alpha_result(k), F_Re_result(k));
            end
        end
    end

    function result = F_CD(gamma,r)    
        result = zeros(size(gamma));

        F_alpha_result = F_alpha(gamma);
        F_Re_result = F_Re(gamma);
        
        for k = 1:length(gamma)
            rad = r(k);     % rのk番目のスカラー値
        
        % 条件分岐
            if rad < 0.066
                result(k) = RootFOIL.CD_alphaRe(F_alpha_result(k), F_Re_result(k));                
            elseif (rad >= 0.066 && rad < 0.073)
                result(k) = ga43.CD_alphaRe(F_alpha_result(k), F_Re_result(k));               
            elseif (rad >= 0.073 && rad < 0.08)
                result(k) = ga42.CD_alphaRe(F_alpha_result(k), F_Re_result(k));                
            elseif (rad >= 0.08 && rad < 0.087)
                result(k) = ga41.CD_alphaRe(F_alpha_result(k), F_Re_result(k));                
            elseif rad >= 0.087 && rad < 0.094
                result(k) = ga40.CD_alphaRe(F_alpha_result(k), F_Re_result(k));                
            elseif rad >= 0.094 && rad < 0.101
                result(k) = ga39.CD_alphaRe(F_alpha_result(k), F_Re_result(k));               
            elseif rad >= 0.101 && rad < 0.108
                result(k) = ga38.CD_alphaRe(F_alpha_result(k), F_Re_result(k));              
            elseif rad >= 0.108 && rad < 0.115
                result(k) = ga37.CD_alphaRe(F_alpha_result(k), F_Re_result(k));                
            elseif rad >= 0.115 && rad < 0.122
                result(k) = ga36.CD_alphaRe(F_alpha_result(k), F_Re_result(k));               
            elseif rad >= 0.122 && rad < 0.129
                result(k) = ga35.CD_alphaRe(F_alpha_result(k), F_Re_result(k));               
            elseif rad >= 0.129 && rad < 0.136
                result(k) = ga34.CD_alphaRe(F_alpha_result(k), F_Re_result(k));             
            elseif rad >= 0.136 && rad < 0.143
                result(k) = ga33.CD_alphaRe(F_alpha_result(k), F_Re_result(k));                
            elseif rad >= 0.143 && rad < 0.15
                result(k) = ga32.CD_alphaRe(F_alpha_result(k), F_Re_result(k));                
            elseif rad >= 0.15 && rad < 0.157
                result(k) = ga31.CD_alphaRe(F_alpha_result(k), F_Re_result(k));               
            elseif rad >= 0.157 && rad < 0.164
                result(k) = ga30.CD_alphaRe(F_alpha_result(k), F_Re_result(k));               
            elseif rad >= 0.164 && rad < 0.172
                result(k) = ga29.CD_alphaRe(F_alpha_result(k), F_Re_result(k));              
            elseif rad >= 0.172 && rad < 0.179
                result(k) = ga28.CD_alphaRe(F_alpha_result(k), F_Re_result(k));              
            elseif rad >= 0.179 && rad < 0.186
                result(k) = ga27.CD_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.186 && rad < 0.193
                result(k) = ga26.CD_alphaRe(F_alpha_result(k), F_Re_result(k));                
            elseif rad >= 0.193 && rad < 0.2
                result(k) = ga25.CD_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.2 && rad < 0.208
                result(k) = ga24.CD_alphaRe(F_alpha_result(k), F_Re_result(k));             
            elseif rad >= 0.208 && rad < 0.215
                result(k) = ga23.CD_alphaRe(F_alpha_result(k), F_Re_result(k));                
            elseif rad >= 0.215 && rad < 0.222
                result(k) = ga22.CD_alphaRe(F_alpha_result(k), F_Re_result(k));               
            elseif rad >= 0.222 && rad < 0.229
                result(k) = ga21.CD_alphaRe(F_alpha_result(k), F_Re_result(k));                
            elseif rad >= 0.229 && rad < 0.236
                result(k) = ga20.CD_alphaRe(F_alpha_result(k), F_Re_result(k));        
            elseif rad >= 0.236 && rad < 0.244
                result(k) = ga19.CD_alphaRe(F_alpha_result(k), F_Re_result(k));               
            elseif rad >= 0.244 && rad < 0.251
                result(k) = ga18.CD_alphaRe(F_alpha_result(k), F_Re_result(k));               
            elseif rad >= 0.251 && rad < 0.258
                result(k) = ga17.CD_alphaRe(F_alpha_result(k), F_Re_result(k));              
            elseif rad >= 0.258 && rad < 0.265
                result(k) = ga16.CD_alphaRe(F_alpha_result(k), F_Re_result(k));              
            elseif rad >= 0.265 && rad < 0.272
                result(k) = ga15.CD_alphaRe(F_alpha_result(k), F_Re_result(k));               
            elseif rad >= 0.272 && rad < 0.279
                result(k) = ga14.CD_alphaRe(F_alpha_result(k), F_Re_result(k));           
            elseif rad >= 0.279 && rad < 0.286
                result(k) = ga13.CD_alphaRe(F_alpha_result(k), F_Re_result(k));           
            elseif rad >= 0.286 && rad < 0.293
                result(k) = ga12.CD_alphaRe(F_alpha_result(k), F_Re_result(k));              
            elseif rad >= 0.293 && rad < 0.3
                result(k) = ga11.CD_alphaRe(F_alpha_result(k), F_Re_result(k));              
            elseif rad >= 0.3 && rad < 0.307
                result(k) = ga10.CD_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.307 && rad < 0.314
                result(k) = ga9.CD_alphaRe(F_alpha_result(k), F_Re_result(k));               
            elseif rad >= 0.314 && rad < 0.321
                result(k) = ga8.CD_alphaRe(F_alpha_result(k), F_Re_result(k));               
            elseif rad >= 0.321 && rad < 0.329
                result(k) = ga7.CD_alphaRe(F_alpha_result(k), F_Re_result(k));                
            elseif rad >= 0.329 && rad < 0.336
                result(k) = ga6.CD_alphaRe(F_alpha_result(k), F_Re_result(k));               
            elseif rad >= 0.336 && rad < 0.343
                result(k) = ga5.CD_alphaRe(F_alpha_result(k), F_Re_result(k));              
            elseif rad >= 0.343 && rad < 0.35
                result(k) = ga4.CD_alphaRe(F_alpha_result(k), F_Re_result(k));               
            elseif rad >= 0.35 && rad < 0.357
                result(k) = ga3.CD_alphaRe(F_alpha_result(k), F_Re_result(k));      
            elseif rad >= 0.357 && rad < 0.364
                result(k) = ga2.CD_alphaRe(F_alpha_result(k), F_Re_result(k));      
            elseif rad >= 0.364 && rad < 0.371
                result(k) = ga1.CD_alphaRe(F_alpha_result(k), F_Re_result(k));
            elseif rad >= 0.371
                result(k) = DAE51.CD_alphaRe(F_alpha_result(k), F_Re_result(k));
            end
        end
    end
    %}

    % caution! 単一翼型か複数翼型でF_CLとF_CDの入力変数が異なるので注意すること。
    F_L = @(gamma) 0.5 * delta_B * RHO * chord .* F_CL(gamma) .* (F_V(gamma).^2);
    F_D = @(gamma) 0.5 * delta_B * RHO * chord .* F_CD(gamma) .* (F_V(gamma).^2);
    %F_L = @(gamma) 0.5 * delta_B * RHO * chord .* F_CL(gamma, r) .* (F_V(gamma).^2);
    %F_D = @(gamma) 0.5 * delta_B * RHO * chord .* F_CD(gamma, r) .* (F_V(gamma).^2);

    F_gamma = @(gamma) F_L(gamma) ./ (RHO .* F_V(gamma) .* delta_B);

    F_T = @(gamma) B * sum(F_L(gamma) .* cos(F_phi(gamma)) - F_D(gamma) .* sin(F_phi(gamma)));
    F_P = @(gamma) B * Omega * sum((F_D(gamma) .* cos(F_phi(gamma)) + F_L(gamma) .* sin(F_phi(gamma))) .* r);
    F_eta = @(gamma) U_inf * F_T(gamma) / F_P(gamma);

    
    % 初期値の計算
    fV = sqrt(U_inf.^2 + (Omega .* r).^2);
    fphi = atan(U_inf ./ (Omega .* r));
    falpha = theta - rad2deg(fphi);
    fRe = fV .* chord / NU;
    fCL = DAE51.CL_alphaRe(falpha, fRe);
    fCD = DAE51.CD_alphaRe(falpha, fRe);
    fL = 0.5 * delta_B * RHO * chord .* fCL .* (fV.^2);
    fD = 0.5 * delta_B * RHO * chord .* fCD .* (fV.^2);
    fT = B * sum(fL .* cos(fphi) - fD .* sin(fphi));
    fP = B * Omega * sum((fD .* cos(fphi) + fL .* sin(fphi)) .* r);
    feta = U_inf * fT / fP;
    fgamma = fL ./ (RHO .* fV .* delta_B);
    

    % 目的関数
    Object = @(gamma) norm(F_gamma(gamma) - gamma);

    % 初期値候補
    f2gamma = sin(linspace(0,1,N)*pi) + 0.01;
    f1gamma = fgamma;

    %二次差分(微分)制約...振動対策
    function [c, ceq] = smooth_penalty(gamma, eps, N)

        e = ones(N,1);
        D2 = spdiags([e -2*e e], -1:1, N, N);
        diff2 = D2*gamma;
        diff2 = diff2(2:end-1);

        % 不等式制約: |diff2| <= eps
        c = [diff2 - eps;   % diff2 <= eps
            -diff2 - eps]; % -diff2 <= eps

        % 等式制約なし
        ceq = [];
    end

    nonlcon = @(x) smooth_penalty(x, eps, N);

    % 初期値選択＆最適化
    options = optimoptions('fmincon', 'MaxIterations',1e8, 'MaxFunctionEvaluations',1e8, 'OptimalityTolerance',1e-30, 'Display', 'iter');
    
    disp('')
    disp('Start exploring')

    try
        % まずf1gammaで試す
        [gamma_opt, ~, exitflag] = fmincon(Object, f1gamma, [], [], [], [], [], [], nonlcon, options);
        if exitflag <= 0 || isnan(F_T(gamma_opt))
            error('failed with f1gamma');
        end
    catch
        % f1gamma失敗ならf2gammaで再試行
        try
            [gamma_opt, ~, exitflag] = fmincon(Object, f2gamma, [], [], [], [], [], [], nonlcon, options);
            if exitflag <= 0 || isnan(F_T(gamma_opt))
                error('failed with f2gamma');
            end
        catch
            warning('解析失敗');
            T=nan; P=nan; eta=nan; norm_val=nan;
            return;
        end
    end

    % 結果
    T = F_T(gamma_opt);
    P = F_P(gamma_opt);
    eta = F_eta(gamma_opt);
    norm_val = norm(F_gamma(gamma_opt) - gamma_opt);

    fprintf('解析結果 @ U_inf=%.2f, RPM=%.2f\n', U_inf, RPM);
    fprintf('推力: %.4f N\n', T);
    fprintf('パワー: %.4f W\n', P);
    fprintf('効率: %.2f %%\n', eta*100);
    fprintf('Norm: %.6f\n', norm_val);

    %%{
    plot(r, gamma_opt, 'b-', 'LineWidth', 2);
    hold on;                             
    plot(r, F_gamma(gamma_opt), 'r--', 'LineWidth', 2);
    hold off;                            
    legend('gammma','Fgamma');
    xlabel('r');
    ylabel('Γ');
    title(sprintf('Comparison of Γ(@ U_inf=%.1f, RPM=%.0f)', U_inf, RPM));
    saveas(gcf,'Comparison of Γ.png')
    %}

    %%{
    plot(r, F_alpha(gamma_opt), 'b-', 'LineWidth', 2);
    xlabel('r');
    ylabel('alpha')
    title(sprintf('迎角分布 (@ U_inf=%.1f, RPM=%.0f)', U_inf, RPM));
    saveas(gcf,'迎角分布.png')
    %}

end

%ある点での性能が知りたいならこれで良い。
analysis_at(7.9, 143, large_table, delta_B, CP, N, ...
    dae51);

%% ロバスト性解析

%{
U_inf_vals = linspace(7,9,3);
RPM_vals = linspace(120,160,3);

[U_grid,RPM_grid] = meshgrid( U_inf_vals, RPM_vals);

result_T = zeros(size(U_grid));
result_P = zeros(size(U_grid));
result_eta = zeros(size(U_grid));
result_norm = zeros(size(U_grid));

for i = 1:length(U_inf_vals)
    for j = 1:length(RPM_vals)
        
        disp('');
        fprintf('Start analysing @ U_inf = %.2f m/s, RPM = %.2f...\n',U_inf_vals(i),RPM_vals(j));
        disp('');

        [T, P, eta, norm_val] = analysis_at( U_inf_vals(i), RPM_vals(j), large_table, delta_B, CP, N, dae51)
        
        result_T(i,j) = T;
        result_P(i,j) = P;
        result_eta(i,j) = eta;
        result_norm(i,j) = norm_val;

        if (result_T(i,j) <= 0) || (result_P(i,j) <= 0)
            result_eta(i,j) = NaN;
        end

        disp('');
    end
end


% 折れ線グラフの重ね合わせ
figure('Position',[100 100 1200 1000])

% U_infに対するTの変化
subplot(2,2,1)
hold on
for i = 1:length(RPM_vals)
    plot(U_inf_vals, result_T(i,:), '-o', 'LineWidth',1, 'MarkerSize',3, ...
        'DisplayName', sprintf('RPM = %.0f', RPM_vals(i)));
end
xlabel('U_{inf} [m/s]')
ylabel('Thrust T [N]')
title('Thrust vs U_{inf}')
legend('Location','best')
grid on
hold off

% U_infに対するPの変化
subplot(2,2,2)
hold on
for i = 1:length(RPM_vals)
    plot(U_inf_vals, result_P(i,:), '-o', 'LineWidth',1, 'MarkerSize',3, ...
        'DisplayName', sprintf('RPM = %.0f', RPM_vals(i)));
end
xlabel('U_{inf} [m/s]')
ylabel('Power')
title('Power vs U_{inf}')
legend('Location','best')
grid on
hold off

% U_infに対するetaの変化
subplot(2,2,3)
hold on
for i = 1:length(RPM_vals)
    plot(U_inf_vals, result_eta(i,:), '-o', 'LineWidth',1, 'MarkerSize',3, ...
        'DisplayName', sprintf('RPM = %.0f', RPM_vals(i)));
end
xlabel('U_{inf} [m/s]')
ylabel('Efficiency')
title('Efficiency vs U_{inf}')
legend('Location','best')
grid on
hold off

% U_infに対するnormの変化
subplot(2,2,4)
hold on
for i = 1:length(RPM_vals)
    plot(U_inf_vals, result_norm(i,:), '-o', 'LineWidth',1, 'MarkerSize',3, ...
        'DisplayName', sprintf('RPM = %.0f', RPM_vals(i)));
end
xlabel('U_{inf} [m/s]')
ylabel('norm')
title('norm vs U_{inf}')
legend('Location','best')
grid on
hold off

saveas(gcf,'Robustness analysis(2D).png')


% 等位線プロット
figure('Position',[100 100 1200 1000])

% 推力分布
subplot(2,2,1)
contourf(U_grid, RPM_grid, result_T, 50, 'LineStyle','none')
colorbar
xlabel('U_{inf}')
ylabel('RPM')
title('Contour Plot of Thrust')

% パワー分布
subplot(2,2,2)
contourf(U_grid, RPM_grid, result_P, 50, 'LineStyle','none')
colorbar
xlabel('U_{inf}')
ylabel('RPM')
title('Contour Plot of Power')

% 効率分布
subplot(2,2,3)
contourf(U_grid, RPM_grid, result_eta, 50, 'LineStyle','none')
colorbar
xlabel('U_{inf}')
ylabel('RPM')
title('Contour Plot of Efficiency')

% norm分布
subplot(2,2,4)
contourf(U_grid, RPM_grid, result_norm, 50, 'LineStyle','none')
colorbar
xlabel('U_{inf}')
ylabel('RPM')
title('Contour Plot of norm')

saveas(gcf,'Robustness analysis(3D).png')
%}

%{
%速度固定ver
U = 7.9;
RPM_vals = linspace(0,200,21);

result_T = zeros(size(RPM_vals));
result_P = zeros(size(RPM_vals));
result_eta = zeros(size(RPM_vals));
result_norm = zeros(size(RPM_vals));

for i =1:length(RPM_vals)
    
    disp('');
    fprintf('Start analysing @ U_inf = %.2f m/s, RPM = %.2f...\n', U, RPM_vals(i));
    disp('');

    [T, P, eta, norm_val] = analysis_at( U, RPM_vals(i), large_table, delta_B, CP, N, dae51)

    result_T(i) = T;
    result_P(i) = P;
    result_eta(i) = eta;
    result_norm(i) = norm_val;

    if (result_T(i) <= 0) || (result_P(i) <= 0)
       result_eta(i) = NaN;
    end

    disp('');
end


figure('Position',[100 100 1200 1000])

% RPMに対するTの変化
subplot(2,2,1)
hold on
plot(RPM_vals, result_T, '-o', 'LineWidth',1, 'MarkerSize',3);
xlabel('RPM')
ylabel('Thrust T [N]')
title('Thrust vs RPM')
grid on
hold off

% RPMに対するPの変化
subplot(2,2,2)
hold on
plot(RPM_vals, result_P, '-o', 'LineWidth',1, 'MarkerSize',3);
xlabel('RPM')
ylabel('Power P [W]')
title('Power vs RPM')
grid on
hold off

% RPMに対するetaの変化
subplot(2,2,3)
hold on
plot(RPM_vals, result_eta, '-o', 'LineWidth',1, 'MarkerSize',3);
xlabel('RPM')
ylabel('eta [%]')
title('eta vs RPM')
grid on
hold off

% U_infに対するnormの変化
subplot(2,2,4)
hold on
plot(RPM_vals, result_norm, '-o', 'LineWidth',1, 'MarkerSize',3);
xlabel('RPM')
ylabel('Norm')
title('Norm vs RPM')
grid on
hold off

saveas(gcf,'Robustness analysis(2D).png')
%}


%回転数固定ver
RPM = 143;
U_inf_vals = linspace(0,12,13);

result_T = zeros(size(U_inf_vals));
result_P = zeros(size(U_inf_vals));
result_eta = zeros(size(U_inf_vals));
result_norm = zeros(size(U_inf_vals));

for i =1:length(U_inf_vals)
    
    disp('');
    fprintf('Start analysing @ U_inf = %.2f m/s, RPM = %.2f...\n', U_inf_vals(i), RPM);
    disp('');

    [T, P, eta, norm_val] = analysis_at( U_inf_vals(i), RPM, large_table, delta_B, CP, N, dae51)

    result_T(i) = T;
    result_P(i) = P;
    result_eta(i) = eta;
    result_norm(i) = norm_val;

    if (result_T(i) <= 0) || (result_P(i) <= 0)
       result_eta(i) = NaN;
    end

    disp('');
end


figure('Position',[100 100 1200 1000])

% U_infに対するTの変化
subplot(2,2,1)
hold on
plot(U_inf_vals, result_T, '-o', 'LineWidth',1, 'MarkerSize',3);
xlabel('U_inf')
ylabel('Thrust T [N]')
title('Thrust vs U_inf')
grid on
hold off

% U_infに対するPの変化
subplot(2,2,2)
hold on
plot(U_inf_vals, result_P, '-o', 'LineWidth',1, 'MarkerSize',3);
xlabel('U_inf')
ylabel('Power P [W]')
title('Power vs U_inf')
grid on
hold off

% U_infに対するetaの変化
subplot(2,2,3)
hold on
plot(U_inf_vals, result_eta, '-o', 'LineWidth',1, 'MarkerSize',3);
xlabel('U_inf')
ylabel('eta [%]')
title('eta vs U_inf')
grid on
hold off

% U_infに対するnormの変化
subplot(2,2,4)
hold on
plot(U_inf_vals, result_norm, '-o', 'LineWidth',1, 'MarkerSize',3);
xlabel('U_inf')
ylabel('Norm')
title('Norm vs U_inf')
grid on
hold off

saveas(gcf,'Robustness analysis(2D).png')
