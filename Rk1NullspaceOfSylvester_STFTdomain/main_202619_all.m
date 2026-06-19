%% 1. 設定と初期化
clear; clc; close all;

fs = 11025
% DGT tool の初期化
F = DGTtool("windowLength", 512, "windowShift", 128);

% List of IRs for evaluation
ir_paths = { ...
    %"logic03sWoodenBooth8000.wav", ...
    "../data/wavIR/IR_sweep_15s_45Hzto22kHz_FS16kHz.v370.wav", ... 
    "../data/wavIR/01.0s Villa Bathroom-OST.wav", ... 
    "../data/BF TL SPACE LIBRARY/Drumbrella/Drumbrella 5'.R.wav"...
};
num_ir = length(ir_paths);

% locating sound file
instrument_dir = "input/Goldberg_clav";
fpaths = dir(fullfile(instrument_dir, '*.wav')); 
num_files = length(fpaths);

if num_files < 2
    error('ペアを作るために、フォルダ内に少なくとも2つのWAVファイルが必要です。');
end

% 全ファイルの組み合わせ（ペア）のインデックスを生成
pair_indices = nchoosek(1:num_files, 2);
num_pairs = size(pair_indices, 1);

fprintf('実験開始: IR数 = %d, 総ファイル数 = %d (全 %d ペア)\n', num_ir, num_files, num_pairs);

% 結果蓄積用のセル配列 (行: IR種別, 列: 各手法のSDR改善量)
% 評価値は各ペアの2つの信号の平均改善量（ΔSDR = 処理後SDR - 観測SDR）を格納します
delta_SDR_CR   = zeros(num_pairs, num_ir);
delta_SDR_Lift = zeros(num_pairs, num_ir);
delta_SDR_WPE  = zeros(num_pairs, num_ir);

% 計算時間蓄積用
time_CR   = zeros(num_pairs, num_ir);
time_Lift = zeros(num_pairs, num_ir);
time_WPE  = zeros(num_pairs, num_ir);

% 箱ひげ図用のラベル集計用（実際のRT60を格納）
rt60_labels = cell(1, num_ir);

%% 2. 巨大総当たりループの実行
for ir_idx = 1:num_ir
    % IRの読み込みとサンプリング・カット処理
    current_ir_path = ir_paths{ir_idx};
    [ir_raw, ~] = audioread(current_ir_path);
    
    rt = reverb_time(ir_raw, fs);
    rt60_labels{ir_idx} = sprintf('RT60 = %.1fs', rt);
    fprintf('\n=========================================\n');
    fprintf('評価中: IR %d/%d (%s)\n', ir_idx, num_ir, rt60_labels{ir_idx});
    fprintf('=========================================\n');
    
    ir = ir_raw(1:floor(fs*rt)); 
    ir = ir(:) / norm(ir);
    IR = F(ir);
    Lir = size(IR, 2);
    
    % CR法用のトリミング位置調整 (お使いのコードのパラメータを継承)
    padding = -2;
    Lir_proc = Lir + padding;

    for p_idx = 1:num_pairs
        idx1 = pair_indices(p_idx, 1);
        idx2 = pair_indices(p_idx, 2);
        fprintf(' -> ペア %d/%d: [%s] & [%s]\n', p_idx, num_pairs, fpaths(idx1).name, fpaths(idx2).name);
        
        % ドライ音の読み込みとノイズ付加
        ss = cell(2,1);
        obs = cell(2,1);
        X = cell(2,1);
        l = zeros(1,2);
        
        % 一時的にWavファイルを保存するパス（Lifting法などの関数へ渡す用）
        tmp_speech_paths = cell(2,1);
        
        for i = 1:2
            f_actual_idx = pair_indices(p_idx, i);
            audiofilename = fullfile(instrument_dir, fpaths(f_actual_idx).name);
            
            % ドライ音読み込み（WPEのfs=16k前提のラッパーに合わせるためfsチェックを推奨）
            [temp, fs_audio] = audioread(audiofilename);
            
            % 擬似的なノイズ付加
            noise = randn(size(temp(:,1)));
            noise = 0.1 * noise * norm(temp(:,1)) / norm(noise);
            ss{i} = temp(:,1) + noise;
            
            % 残響の畳み込み
            obs{i} = conv(ss{i}, ir);
            X{i} = F(obs{i});
            
            % 後段のLifting法にクリーン音をパス経由で渡すための一時保存
            tmp_speech_paths{i} = sprintf('tmp_speech_p%d_ch%d.wav', p_idx, i);
            audiowrite(tmp_speech_paths{i}, ss{i}, fs);
        end
        
        Fq = size(X{1}, 1);
        N_len = [size(X{1},2), size(X{2},2)];
        for i = 1:2, l(i) = N_len(i) - Lir_proc + 1; end
        
        %% --- 手法A: 提案（原始的CR法） ---
        RET1 = zeros(Fq, l(1)); RET2 = zeros(Fq, l(2));
        tic
        parfor f = 1:Fq
            X1f = X{1}(f, :)';
            X2f = X{2}(f, :)';
            Xf = [X1f; X2f];
            XfNorm = norm(Xf);
            
            Sylvester = [convmtx(X2f, l(1)) convmtx(-X1f, l(2))];
            [~, ~, V] = svd(Sylvester, "econ");
            Vf = V(:,end);
            Vf = XfNorm * Vf / norm(Vf);
            
            RET1(f,:) = Vf(1:l(1))';
            RET2(f,:) = Vf(l(1)+1:end)';
        end
        ret_cr = cell(2,1);
        ret_cr{1} = F.pinv(RET1);  
        ret_cr{2} = F.pinv(RET2);
        time_CR(p_idx, ir_idx) = toc;
        
        %% --- 手法B: Lifting法 ---
        ret_lift = cell(2,1);
        tic
        for i = 1:2
            RET_lift_sub = dereverb(tmp_speech_paths{i}, current_ir_path);
            ret_lift{i} = F.pinv(RET_lift_sub);
        end
        time_Lift(p_idx, ir_idx) = toc;
        
        %% --- 手法C: WPE ---
        ret_wpe = cell(2,1);
        tic
        for i = 1:2
            ret_wpe{i} = WPE1c_wrapper(obs{i}, fs);
        end
        time_WPE(p_idx, ir_idx) = toc;
        
        %% --- 評価（SDR算出） ---
        pair_sdr_obs  = zeros(2,1);
        pair_sdr_cr   = zeros(2,1);
        pair_sdr_lift = zeros(2,1);
        pair_sdr_wpe  = zeros(2,1);
        
        for i = 1:2
            % 長さのミスマッチを最小に揃える
            len_ss = length(ss{i});
            len_cr = length(ret_cr{i});
            len_lift = length(ret_lift{i});
            len_wpe = length(ret_wpe{i});
            
            min_len = min([len_ss, len_cr, len_lift, len_wpe]);
            
            % 各信号のクリッピング
            s_target = ss{i}(1:min_len)';
            s_obs    = obs{i}(1:min_len)';
            s_cr     = ret_cr{i}(1:min_len)';
            s_lift   = ret_lift{i}(1:min_len)';
            s_wpe    = ret_wpe{i}(1:min_len)';
            
            % SDR計算
            pair_sdr_obs(i)  = bss_eval_sources(s_obs,  s_target);
            pair_sdr_cr(i)   = bss_eval_sources(s_cr,   s_target);
            pair_sdr_lift(i) = bss_eval_sources(s_lift, s_target);
            pair_sdr_wpe(i)  = bss_eval_sources(s_wpe,  s_target);
        end
        
        % ペア内の2チャネルの平均SDR改善量 (ΔSDR) を蓄積
        delta_SDR_CR(p_idx, ir_idx)   = mean(pair_sdr_cr - pair_sdr_obs);
        delta_SDR_Lift(p_idx, ir_idx) = mean(pair_sdr_lift - pair_sdr_obs);
        delta_SDR_WPE(p_idx, ir_idx)  = mean(pair_sdr_wpe - pair_sdr_obs);
        
        % 一時ファイルの削除（ディスク逼迫を防止）
        for i = 1:2
            if exist(tmp_speech_paths{i}, 'file'), delete(tmp_speech_paths{i}); end
        end
    end
end

%% 3. 箱ひげ図 (Boxchart) の自動プロット（重複エラー回避版）
fprintf('\n実験完了。プロットを作成中...\n');

% boxchart用にデータをテーブル形式に整理
total_samples = num_pairs * num_ir * 3; % 3手法分
plot_SDR = zeros(total_samples, 1);
plot_IR  = cell(total_samples, 1);
plot_Method = cell(total_samples, 1);

% 重複を確実に防ぐため、ValueSet（マスターのラベルリスト）を一意に生成
unique_labels = cell(1, num_ir);
for ir_idx = 1:num_ir
    % ラベルの先頭に「IR 1: 」のように番号を振ることで、文字の重複を完全に防ぎます
    unique_labels{ir_idx} = sprintf('IR %d (%s)', ir_idx, rt60_labels{ir_idx});
end

idx = 1;
for ir_idx = 1:num_ir
    for p_idx = 1:num_pairs
        % CR法
        plot_SDR(idx) = delta_SDR_CR(p_idx, ir_idx);
        plot_IR{idx}  = unique_labels{ir_idx};
        plot_Method{idx} = 'Proposed CR';
        idx = idx + 1;
        
        % Lifting法
        plot_SDR(idx) = delta_SDR_Lift(p_idx, ir_idx);
        plot_IR{idx}  = unique_labels{ir_idx};
        plot_Method{idx} = 'Lifting';
        idx = idx + 1;
        
        % WPE法
        plot_SDR(idx) = delta_SDR_WPE(p_idx, ir_idx);
        plot_IR{idx}  = unique_labels{ir_idx};
        plot_Method{idx} = 'Monaural WPE';
        idx = idx + 1;
    end
end

% カテゴリカル変数に変換して表示順序を固定
plot_IR = categorical(plot_IR, unique_labels);
plot_Method = categorical(plot_Method, {'Proposed CR', 'Lifting', 'Monaural WPE'});

% plot Graph
figure('Position', [100, 100, 750, 480]);
h = boxchart(plot_IR, plot_SDR, 'GroupByColor', plot_Method);

% グラフの調整
grid on;
ylabel('\Delta SDR [dB]', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Room Impulse Response (RIR)', 'FontSize', 12);
title('Dereverberation Performance Across Different RT60', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 10);
set(gca, 'FontSize', 11);

% 平均処理時間の表示
fprintf('\n---- mean Ex time (per pair) ----\n')
fprintf('Proposed CR : %.4f 秒\n', mean(time_CR(:)));
fprintf('Lifting     : %.4f 秒\n', mean(time_Lift(:)));
fprintf('Monaural WPE: %.4f 秒\n', mean(time_WPE(:)));