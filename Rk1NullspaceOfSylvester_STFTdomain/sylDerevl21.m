function [v1, v2, info] = sylDerevl21(X1f, X2f, l, rho, alpha)
Xf = [X1f; X2f];
XfNorm = norm(Xf);   % XfPhase = Xf ./ abs(Xf);

% Null-space calculation via (not structured) low-rank approximation
Sylvester = [convmtx(X2f, l(1)) convmtx(-X1f, l(2))];
[~, ~, V] = svd(Sylvester, "econ");

% 事前準備
AtA = Sylvester' * Sylvester;
mu = 10^-6;
decompAtA = decomposition(AtA + (rho + mu) * eye(size(AtA)));

% 変数の初期化
u_init = V(:, end);
u = u_init;
z = u_init;
y = zeros(size(u));

% --- Definition of Group ---
group_size = 4;
num_elements = length(u);
num_groups = ceil(num_elements / group_size);

for iter = 1:20
    % 1. u-update
    u = decompAtA \ (rho * z - y + mu * u_init);

    % 2. z-update (Block Soft-Thresholding)
    v = u + y / rho;
    z = zeros(size(v)); % 詰め替え用のバッファ

    for g = 1:num_groups
        % 各グループのインデックスを抽出
        idx = ((g-1)*group_size + 1) : min(g*group_size, num_elements);
        v_g = v(idx);

        % グループの L2 ノルムを計算
        norm_vg = norm(v_g, 2);

        % Block Soft-Thresholding
        if norm_vg > (alpha / rho)
            % グループ内の全要素を、共通の縮小係数でスケールダウン
            z(idx) = v_g * (1 - (alpha / rho) / norm_vg);
        else
            % グループのエネルギーが閾値以下なら、グループ丸ごとゼロに
            z(idx) = 0;
        end
    end

    if norm(z) < 1e-12
        z = u_init;
    else
        z = z / norm(z);
    end

    % 3. LargangeMultiplier-update
    y = y + rho * (u - z);
end

% 最終的なクリーン信号の推定値として z を使用する
Vf = z;    Vf = XfNorm * Vf / norm(Vf); % scale recovery

% % subband IR estimation
% A = [convmtx(Vf(1:l(1)), L(1)); convmtx(Vf(l(1)+1:end), L(1))];     b = [X1f; X2f];
% Hhat = A \ b; % calc. IR (TF domain)
% HhatPhase = Hhat / abs(Hhat);
% [Hhatmax, idx_HhatPeak] = max(abs(Hhat));

v1 = Vf(1:l(1));   v2 = Vf(l(1)+1:end);
end