function d = efficient_deconvolution(y, h, L, Q, options)
    P = length(h);
    L_P = P + Q;
    L_a = L + L_P;
    M = length(y);

    % Step 1: Precompute H_L
    H_L = fft(h, L);

    % 正則化と閾値処理
    if options.adaptiveThreshold
        epsilon = 0.01 * max(abs(H_L));  % 動的設定
    else
        epsilon = options.epsilon;       % 固定値
    end

    if options.regularization
        % 高度な正則化（Tikhonov型）
        H_L_safe = conj(H_L) ./ (abs(H_L).^2 + epsilon);
    else
        % 単純なしきい値処理
        H_L_safe = H_L;
        H_L_safe(abs(H_L) < epsilon) = epsilon;
    end

    % ブロック処理
    num_blocks = floor((M - L) / L_a) + 1;
    d = zeros(1, num_blocks * L);

    for r = 0:num_blocks-1
        idx_start = r * L_a + 1;
        idx_end = idx_start + L - 1;
        if idx_end > M
            break;
        end

        y_block = y(idx_start:idx_end);
        y_rL = y_block .* ones(L, 1);  % 方形窓

        Y_rL = fft(y_rL, L);

        if options.regularization
            D_r = Y_rL .* H_L_safe;
        else
            D_r = Y_rL ./ H_L_safe;
        end

        d_r = real(ifft(D_r, L));
        d_idx = r * L + (1:L);
        d(d_idx) = d(d_idx) + d_r';
    end
end

% options = struct( ...
%     'regularization', true, ...       % 高度な正則化を使うかどうか
%     'adaptiveThreshold', true, ...    % 閾値を動的に設定するかどうか
%     'epsilon', 1e-3 ...               % 固定閾値（adaptiveThreshold=falseのとき使用）
% );
