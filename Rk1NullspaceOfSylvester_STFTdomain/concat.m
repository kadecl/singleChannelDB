% 擬似的に位相をバカにさせた推定IRを作る
H_true = fft(h_true);
random_phase = exp(j * 2 * pi * rand(size(H_true))); % ランダム位相
H_distorted = abs(H_true) .* random_phase;          % 振幅は無傷、位相は崩壊
h_distorted = real(ifft(H_distorted));               % グチャグチャな時間IR

% これを実数ケプストラムに変換（位相ノイズが消え、対称な形状が得られる）
c_input = rceps(h_distorted);