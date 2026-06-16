%% ========================================================================
%% 1D近接作用素：L2,1模倣学習データ生成 ＆ 再トレーニング
%% ========================================================================
T_total = l(1) + l(2); 
num_samples = 5000; 

X_train = zeros(T_total, 2, 1, num_samples, 'single'); % 入力: v
Y_train = zeros(T_total, 2, 1, num_samples, 'single'); % 教師: L2,1を通した z

% グループスパースのパラメータ（吉野さんの元の設定）
group_size = 4;
num_elements = T_total;
num_groups = ceil(num_elements / group_size);

for i = 1:num_samples
    % ADMMループ内で発生するリアルな暫定ベクトル v をシミュレート
    % (適度な残響感とノイズを混ぜる)
    v_raw = randn(T_total, 1) + 1j * randn(T_total, 1);
    v_raw = filter(0.8.^(0:10), 1, v_raw); % 残響ブラー
    v_raw = v_raw / max(abs(v_raw)); % 最大値を1に正規化
    
    % --- 従来の吉野さんのグループスパース（Block Soft-Thresholding）をそのまま実行 ---
    z_true = zeros(size(v_raw));
    for g = 1:num_groups
        idx = ((g-1)*group_size + 1) : min(g*group_size, num_elements);
        v_g = v_raw(idx);
        norm_vg = norm(v_g, 2);
        
        % alpha/rho の閾値。入力が1に正規化されているので、少し強めに設定（例: 0.15）
        th = 0.15; 
        if norm_vg > th
            z_true(idx) = v_g * (1 - th / norm_vg);
        else
            z_true(idx) = 0;
        end
    end
    
    % 入力（グループスパースをかける前）
    X_train(:, 1, 1, i) = real(v_raw);
    X_train(:, 2, 1, i) = imag(v_raw);
    
    % 教師（グループスパースをかけた後の、美しいスパース波形）
    Y_train(:, 1, 1, i) = real(z_true);
    Y_train(:, 2, 1, i) = imag(z_true);
end

%% ========================================================================
%% 2. 1D-CNN ネットワーク構築 (修正版)
%% ========================================================================
% 入力サイズを [T_total 2 1] に明示的に固定します
%% 2. 1D-ResNet ネットワーク構築 (強力版)
% 入力: [T_total x 2 x 1]
lgraph = layerGraph();

% メインパスの定義
mainLayers = [
    imageInputLayer([T_total 2 1], 'Name', 'input', 'Normalization', 'none')
    
    % 第1ブロック
    convolution2dLayer([7 1], 64, 'Padding', 'same', 'Name', 'conv1')
    reluLayer('Name', 'relu1')
    
    % 第2ブロック (ここから残差構造を意識)
    convolution2dLayer([5 1], 64, 'Padding', 'same', 'Name', 'conv2')
    reluLayer('Name', 'relu2')
    convolution2dLayer([5 1], 64, 'Padding', 'same', 'Name', 'conv3')
    
    % 加算層（スキップ接続用）
    additionLayer(2, 'Name', 'add1')
    reluLayer('Name', 'relu3')
    
    % 第3ブロック
    convolution2dLayer([3 1], 64, 'Padding', 'same', 'Name', 'conv4')
    reluLayer('Name', 'relu4')
    convolution2dLayer([3 1], 64, 'Padding', 'same', 'Name', 'conv5')
    
    additionLayer(2, 'Name', 'add2')
    reluLayer('Name', 'relu5')
    
    % 出力層：1x1畳み込みで [T_total x 2 x 1] にブレンド
    convolution2dLayer([1 1], 1, 'Padding', 'same', 'Name', 'conv_out')
    regressionLayer('Name', 'output')
];

lgraph = addLayers(lgraph, mainLayers);

% スキップ接続（Identity Mapping）の追加
% 入力に近い特徴を、後ろの層に直接足し合わせることで、勾配消失を防ぎ深く学習させます
lgraph = connectLayers(lgraph, 'relu1', 'add1/in2');
lgraph = connectLayers(lgraph, 'relu3', 'add2/in2');

% 学習オプション
options = trainingOptions('adam', ...
    'MaxEpochs', 60, ...
    'MiniBatchSize', 128, ...
    'InitialLearnRate', 1e-3, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 15, ...
    'LearnRateDropFactor', 0.5, ...
    'Shuffle', 'every-epoch', ...
    'Plots', 'training-progress', ...
    'ExecutionEnvironment', 'gpu');

net_1d_prox_2 = trainNetwork(X_train, Y_train, lgraph, options);
save('net_1d_prox.mat', 'net_1d_prox_2');