%% ========================================================================
%% 1D近接作用素：L2,1模倣学習データ生成 ＆ 再トレーニング
%% ========================================================================
%% ========================================================================
%% 1D振幅マスキング（Learned Proximal Mask）データ生成 ＆ 訓練
%% ========================================================================
T_total = l(1) + l(2); 
num_samples = 100000; % 30分間フルパワー用に10万件

% 入力：[T_total x 1 x 1 x num_samples]（振幅のみなのでチャンネル数は 1）
X_train = zeros(T_total, 1, 1, num_samples, 'single'); % 入力: 滲んだ振幅
Y_train = zeros(T_total, 1, 1, num_samples, 'single'); % 教師: 理想の1Dマスク(0~1)

for i = 1:num_samples
    % 音響ライクなスパース複素信号の生成
    s_clean_real = randn(T_total, 1) .* (rand(T_total, 1) > 0.85) .* filter(1, [1 -0.9], randn(T_total, 1)); 
    s_clean_imag = randn(T_total, 1) .* (rand(T_total, 1) > 0.85) .* filter(1, [1 -0.9], randn(T_total, 1));
    s_clean = s_clean_real + 1j * s_clean_imag;
    
    % ランダムな残響テイルの付与
    current_decay = 0.70 + rand()*0.22;
    kernel_time = current_decay.^(0:20); kernel_time = kernel_time / sum(kernel_time);
    s_blurred = filter(kernel_time, 1, s_clean);
    
    % ランダムなADMMノイズの付与
    noise = (randn(T_total, 1) + 1j*randn(T_total, 1)) * (0.01 + rand()*0.09);
    s_blurred_noisy = s_blurred + noise;
    
    % --- 【ここが核心】振幅ドメインへの変換とマスクの計算 ---
    mag_clean = abs(s_clean);
    mag_blurred = abs(s_blurred_noisy);
    
    % 入力（滲んだ振幅）を最大値1に正規化
    v_max = max(mag_blurred); if v_max < 1e-12, v_max = 1; end
    mag_blurred_norm = mag_blurred / v_max;
    
    % 理想マスクの計算（クリーン振幅 / 滲み振幅）
    % 0~1の範囲にクリッピングして、完全に安定した教師信号を作ります
    mask_true = mag_clean ./ (mag_blurred + eps);
    mask_true = min(max(mask_true, 0), 1); 
    
    % 格納（3次元目を1、4次元目をサンプルインデックスに）
    X_train(:, 1, 1, i) = single(mag_blurred_norm);
    Y_train(:, 1, 1, i) = single(mask_true);
end

%================================================================================================================
%% 2. 1Dマスク予測用 DeepResNet 構築
lgraph = layerGraph();
num_filters = 96;

mainLayers = [
    imageInputLayer([T_total 1 1], 'Name', 'input', 'Normalization', 'none') % チャンネル数は1
    
    convolution2dLayer([7 1], num_filters, 'Padding', 'same', 'Name', 'conv_init')
    reluLayer('Name', 'relu_init')
    
    % ResBlock 1
    convolution2dLayer([5 1], num_filters, 'Padding', 'same', 'Name', 'res1_c1')
    reluLayer('Name', 'relu1_c1')
    convolution2dLayer([5 1], num_filters, 'Padding', 'same', 'Name', 'res1_c2')
    additionLayer(2, 'Name', 'add1')
    reluLayer('Name', 'relu_res1_out')
    
    % ResBlock 2
    convolution2dLayer([5 1], num_filters, 'Padding', 'same', 'Name', 'res2_c1')
    reluLayer('Name', 'relu2_c1')
    convolution2dLayer([5 1], num_filters, 'Padding', 'same', 'Name', 'res2_c2')
    additionLayer(2, 'Name', 'add2')
    reluLayer('Name', 'relu_res2_out')
    
    % 1x1畳み込みで1チャンネルに集約
    convolution2dLayer([1 1], 1, 'Padding', 'same', 'Name', 'conv_out')
    
    % 【超重要】出力を強制的に 0 〜 1 のマスクの範囲に縛り付ける
    % これにより、無音区間で負の値や異常な高値が出るのを物理的に防ぎます
    sigmoidLayer('Name', 'sigmoid_out') 
    
    regressionLayer('Name', 'output')
];

lgraph = addLayers(lgraph, mainLayers);
lgraph = connectLayers(lgraph, 'relu_init',     'add1/in2');
lgraph = connectLayers(lgraph, 'relu_res1_out', 'add2/in2');

% 30分フルパワー訓練オプション（前回の設定を維持）
options = trainingOptions('adam', ...
    'MaxEpochs', 100, ...
    'MiniBatchSize', 256, ...
    'InitialLearnRate', 1e-3, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 25, ...
    'LearnRateDropFactor', 0.5, ...
    'Shuffle', 'every-epoch', ...
    'Plots', 'training-progress', ...
    'ExecutionEnvironment', 'gpu');

net_1d_mask = trainNetwork(X_train, Y_train, lgraph, options);
save('net_1d_mask.mat', 'net_1d_mask');