function z = modifiedSpectralSubtraction(x, fs, optSTFT)
if nargin < 3
    Y = stft(x, fs);
else
    Y = stft(x, fs, optSTFT);
end
% Y のサイズを取得
[W, L] = size(Y);
Yabs = abs(Y);

% 初期化
alpha = zeros(W, L);

for i = 1:L
    % 各周波数ビン w に対して計算
    for w = 1:W
        % フレーム方向の有効範囲を設定
        m_range = (i+1):L; % 遅延 i を考慮したフレーム範囲
        % 分子と分母を計算
        numerator = abs(Y(w, m_range) .* conj(Y(w, m_range - i))); 
        % |Y(w, m) * conj(Y(w, m-i))|
        denominator = Yabs(w, m_range - i).^2; % |Y(w, m-i)|^2
        % フレーム方向の平均 (期待値)
        alpha(w,i) = mean(numerator ./ denominator);
    end
end

P = zeros(size(Y));
for w = 1:W
    for m = 1:L
        range_i = 1:m - 1;
        P(w,m) = sum(abs(alpha(w,range_i)).^2 .* Yabs(w,m-range_i).^2);
    end
end

G = sqrt( max(1 - P ./ (Yabs.^2 + eps) ,0) );
z = istft(G .* Y);
end
