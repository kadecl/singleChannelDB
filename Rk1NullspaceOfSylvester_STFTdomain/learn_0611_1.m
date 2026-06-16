%% ========================================================================
%% 1D近接作用素（Learned Proximal）学習スクリプト (エラー修正版)
%% ========================================================================

T_total = l(1) + l(2); 
num_samples = 5000; 

% 【修正点】4次元配列として初期化します
% [時間(T_total) x チャンネル(2) x 1(固定) x サンプル数(num_samples)]
X_train = zeros(T_total, 2, 1, num_samples, 'single'); 
Y_train = zeros(T_total, 2, 1, num_samples, 'single'); 

% 時間方向の残響滲みカーネル
decay_rate = 0.85; 
kernel_time = decay_rate.^(0:15); 
kernel_time = kernel_time / sum(kernel_time);

for i = 1:num_samples
    % クリーンな1次元信号の模擬
    s_clean_real = randn(T_total, 1) .* (sin((1:T_total)'/5) > 0); 
    s_clean_imag = randn(T_total, 1) .* (sin((1:T_total)'/5) > 0);
    s_clean = s_clean_real + 1j*s_clean_imag;
    
    % 残響滲みの付与
    s_blurred = filter(kernel_time, 1, s_clean);
    
    % 複素雑音の模擬
    noise = (randn(T_total, 1) + 1j*randn(T_total, 1)) * 0.05;
    s_blurred_noisy = s_blurred + noise;
    
    % 【修正点】4次元目のインデックス i に対して、3次元目を 1 固定にして格納
    X_train(:, 1, 1, i) = real(s_blurred_noisy);
    X_train(:, 2, 1, i) = imag(s_blurred_noisy);
    
    Y_train(:, 1, 1, i) = real(s_clean);
    Y_train(:, 2, 1, i) = imag(s_clean);
end

%% ========================================================================
%% 2. 1D-CNN ネットワーク構築 (修正版)
%% ========================================================================
% 入力サイズを [T_total 2 1] に明示的に固定します
% 2. 1D-CNN ネットワーク構築 
layers = [
    imageInputLayer([T_total 2 1], 'Name', 'input', 'Normalization', 'none')
    
    % 時間方向（縦方向のみ）に畳み込みを行い、特徴を抽出
    convolution2dLayer([5 1], 32, 'Padding', 'same', 'Name', 'conv1')
    reluLayer('Name', 'relu1')
    
    convolution2dLayer([5 1], 64, 'Padding', 'same', 'Name', 'conv2')
    reluLayer('Name', 'relu2')
    
    convolution2dLayer([3 1], 64, 'Padding', 'same', 'Name', 'conv3')
    reluLayer('Name', 'relu3')
    
    % --- 【ここを修正】 ---
    % チャンネル数を増やすのではなく、最後の特徴マップから 
    % 入力と同じ [T_total x 2 x 1] の形状へジャストフィットさせるための畳み込み
    convolution2dLayer([1 1], 1, 'Padding', 'same', 'Name', 'conv_out')
    
    regressionLayer('Name', 'output')
];

% 学習オプション
options = trainingOptions('adam', ...
    'MaxEpochs', 25, ...
    'MiniBatchSize', 64, ...
    'InitialLearnRate', 1e-3, ...
    'Shuffle', 'every-epoch', ...
    'Plots', 'training-progress', ...
    'ExecutionEnvironment', 'gpu');

net_1d_prox = trainNetwork(X_train, Y_train, layers, options);
save('net_1d_prox.mat', 'net_1d_prox');