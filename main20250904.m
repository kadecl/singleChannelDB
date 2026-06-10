fsRe = 8000;
[s1org, fs] = audioread("Rk1NullspaceOfSylvester_STFTdomain/input2/AGt/A.Gt_1.wav");
if size(s1org, 2) > 1, s1org = sum(s1org, 2); end
s1 = resample(s1org, fsRe, fs);
[s2org, fs] = audioread("Rk1NullspaceOfSylvester_STFTdomain/input2/AGt/A.Gt_2.wav");
if size(s2org, 2) > 1, s2org = sum(s2org, 2); end
s2 = resample(s2org, fsRe, fs);


if numel(s1) > numel(s2), s2 = buffer(s2, numel(s1)); 
else, s1 = buffer(s1, numel(s2)); end
Ls = [length(s1) length(s2)];

% ir generating
duration = 0.4;
t = linspace(0, duration, fsRe * duration);  % 時間軸

% 指数的減衰（単純なモデル）
decay_rate = 5;          % 減衰率（大きいほど早く減衰）
ir = exp(-decay_rate * t);

% 最初にインパルスを加える
ir(1) = 1; ir = ir / norm(ir);
Lh = length(ir);

x1 = conv(s1, ir); x2 = conv(s2, ir);
%%
lambda = 10^2;
[s1_hat, s2_hat, h_hat] = gcdDeconv(x1, x2, Ls, Lh, lambda);