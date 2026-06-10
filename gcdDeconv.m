function [s1_hat, s2_hat, h_hat] = gcdDeconv(x1, x2, Ls, Lh, lambda)
    % x1, x2: 観測信号
    % Ls: s1, s2の長さ
    % Lh: hの長さ
    % lambda: 正則化係数

    % 初期化
    s1_init = randn(Ls(1), 1);
    s2_init = randn(Ls(2), 1);
    h_init  = randn(Lh, 1);

    % 重み関数（エネルギー減衰型）
    w = linspace(1, 10, Lh)';  % 時間とともに増加する重み

    % 目的関数
    cost_func = @(params) objective(params, x1, x2, Ls, Lh, lambda, w);

    % 最適化
    options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton', 'UseParallel',true);
    params_init = [s1_init; s2_init; h_init];
    params_hat = fminunc(cost_func, params_init, options);

    % 結果の分離
    s1_hat = params_hat(1:Ls(1));
    s2_hat = params_hat(Ls(1)+1:sum(Ls));
    h_hat  = params_hat(sum(Ls)+1:end);
end

function J = objective(params, x1, x2, Ls, Lh, lambda, w)
    s1 = params(1:Ls(1));
    s2 = params(Ls(1)+1:sum(Ls));
    h  = params(sum(Ls)+1:end);

    x1_hat = conv(s1, h);
    x2_hat = conv(s2, h);

    % 誤差項（必要に応じて長さを調整）
    err1 = x1(1:length(x1_hat)) - x1_hat;
    err2 = x2(1:length(x2_hat)) - x2_hat;

    % 正則化項
    reg = lambda * sum(w .* (h.^2));

    % 総コスト
    J = sum(err1.^2) + sum(err2.^2) + reg;
end
