% 1. ネットワーク層の定義
L = 2048; % IRのサンプル長
layers = [
    sequenceInputLayer(L, 'Name', 'input')
    
    convolution1dLayer(5, 64, 'Padding', 'same')
    reluLayer
    
    convolution1dLayer(5, 128, 'Padding', 'same')
    reluLayer
    
    convolution1dLayer(5, 64, 'Padding', 'same')
    reluLayer
    
    convolution1dLayer(5, 1, 'Padding', 'same') % 出力は1チャネル（IR波形）
    regressionLayer('Name', 'output')
];

% 2. 学習オプションの設定（RTX 4090を指定）
options = trainingOptions('adam', ...
    'MaxEpochs', 50, ...
    'MiniBatchSize', 64, ...
    'InitialLearnRate', 1e-3, ...
    'Plots', 'training-progress', ...
    'ExecutionEnvironment', 'gpu'); % RTX 4090で実行

% 3. 学習の実行
net = trainNetwork(X_train, Y_train, layers, options);