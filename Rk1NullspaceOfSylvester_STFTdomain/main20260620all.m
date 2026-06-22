%% 1. 設定と初期化
clear; clc; close all;

% DGT tool の初期化
F = DGTtool("windowLength", 512, "windowShift", 128);
fs_glob = 11025;

% 評価対象のインパルス応答（IR）のリスト
ir_paths = dir("./ir");
isnt_ir_file = {ir_paths.isdir};
isnt_ir_file = cell2mat(isnt_ir_file);
ir_paths = "./ir/" + {ir_paths(~isnt_ir_file).name}
num_ir = length(ir_paths);

% 親フォルダの指定とサブフォルダ（1〜10）の定義
base_dir = "input/inst";
sub_dirs = dir(base_dir);
dnames = {sub_dirs.name};
isDot = startsWith(dnames, '.');   % 先頭が '.' のものを検出
sub_dirs = sub_dirs(~isDot);
num_sets = length(sub_dirs);

fprintf('実験開始: IR数 = %d, 総データセット（フォルダ）数 = %d\n', num_ir, num_sets);

% 結果蓄積用の配列 (行: 1〜10のセット, 列: 各IR種別)
% 重複がないため、各セットの2信号の平均改善量（ΔSDR）を1つのクリーンなデータとして格納します
delta_SDR_CR   = zeros(num_sets, num_ir);
delta_SDR_Lift = zeros(num_sets, num_ir);
delta_SDR_WPE  = zeros(num_sets, num_ir);

% 計算時間蓄積用
time_CR   = zeros(num_sets, num_ir);
time_Lift = zeros(num_sets, num_ir);
time_WPE  = zeros(num_sets, num_ir);

% 箱ひげ図用のラベル集計用
rt60_labels = cell(1, num_ir);

%% 2. 処理ループ（IR × サブフォルダ）
for ir_idx = 1:num_ir
    % IRの読み込みとサンプリング・カット処理
    current_ir_path = ir_paths{ir_idx};
    [ir_raw, fs] = audioread(current_ir_path);
    
    rt = reverb_time(ir_raw, fs_glob);
    rt60_labels{ir_idx} = sprintf('RT60 = %.1fs', rt);
    fprintf('\n=========================================\n');
    fprintf('評価中: IR %d/%d (%s)\n', ir_idx, num_ir, rt60_labels{ir_idx});
    fprintf('=========================================\n');
    
    ir = ir_raw(1:floor(fs_glob*rt)); 
    ir = ir(:) / norm(ir);
    IR = F(ir);
    Lir = size(IR, 2);
    
    padding = -2;
    Lir_proc = Lir + padding;

    for s_idx = 1:num_sets
        current_sub_dir = fullfile(base_dir, sub_dirs(s_idx).name);
        fpaths = dir(fullfile(current_sub_dir, '*.wav'));
        
        if length(fpaths) < 2
            warning('フォルダ「%s」内にWAVファイルが2つ未満のためスキップします。', sub_dirs(s_idx).name);
            continue;
        end
        
        fprintf(' -> フォルダ 「%s」を処理中...\n', sub_dirs(s_idx).name);
        
        % ドライ音の読み込みとノイズ付加
        ss = cell(2,1);
        obs = cell(2,1);
        X = cell(2,1);
        l = zeros(1,2);
        tmp_speech_paths = cell(2,1);
        
        for i = 1:2
            audiofilename = fullfile(current_sub_dir, fpaths(i).name);
            [temp, fs_audio] = audioread(audiofilename);
            
            % 擬似的なノイズ付加
            noise = randn(size(temp(:,1)));
            noise = 0.1 * noise * norm(temp(:,1)) / norm(noise);
            ss{i} = temp(:,1) + noise;
            
            % 残響の畳み込み
            obs{i} = conv(ss{i}, ir);
            X{i} = F(obs{i});
            
            % 後段のLifting法へ渡す用の一時ファイル（フォルダ直下に保存して競合回避）
            tmp_speech_paths{i} = fullfile(current_sub_dir, sprintf('tmp_speech_ch%d.wav', i));
            audiowrite(tmp_speech_paths{i}, ss{i}, fs_glob);
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
        time_CR(s_idx, ir_idx) = toc;
        
        %% --- 手法B: Lifting法 ---
        ret_lift = cell(2,1);
        tic
        for i = 1:2
            RET_lift_sub = dereverb(X{i}, Lir);
            ret_lift{i} = F.pinv(RET_lift_sub);
        end
        time_Lift(s_idx, ir_idx) = toc;
        
        %% --- 手法C: WPE ---
        ret_wpe = cell(2,1);
        tic
        for i = 1:2
            ret_wpe{i} = WPE1c_wrapper(obs{i}, fs_glob, F, Lir);
        end
        time_WPE(s_idx, ir_idx) = toc;
        
        %% --- 評価（SDR算出） ---
        pair_sdr_obs  = zeros(2,1);
        pair_sdr_cr   = zeros(2,1);
        pair_sdr_lift = zeros(2,1);
        pair_sdr_wpe  = zeros(2,1);
        
        for i = 1:2
            min_len = min([length(ss{i}), length(ret_cr{i}), length(ret_lift{i}), length(ret_wpe{i})]);
            
            s_target = ss{i}(1:min_len)';
            s_obs    = obs{i}(1:min_len)';
            s_cr     = ret_cr{i}(1:min_len)';
            s_lift   = ret_lift{i}(1:min_len)';
            s_wpe    = ret_wpe{i}(1:min_len)';
            
            pair_sdr_obs(i)  = bss_eval_sources(s_obs,  s_target);
            pair_sdr_cr(i)   = bss_eval_sources(s_cr,   s_target);
            pair_sdr_lift(i) = bss_eval_sources(s_lift, s_target);
            pair_sdr_wpe(i)  = bss_eval_sources(s_wpe,  s_target);
        end
        
        % セット内の2信号の平均SDR改善量 (ΔSDR) を蓄積
        delta_SDR_CR(s_idx, ir_idx)   = mean(pair_sdr_cr - pair_sdr_obs);
        delta_SDR_Lift(s_idx, ir_idx) = mean(pair_sdr_lift - pair_sdr_obs);
        delta_SDR_WPE(s_idx, ir_idx)  = mean(pair_sdr_wpe - pair_sdr_obs);
        
        % 一時ファイルの削除
        for i = 1:2
            if exist(tmp_speech_paths{i}, 'file'), delete(tmp_speech_paths{i}); end
        end
    end
end

%% 3. 箱ひげ図 (Boxchart) の自動プロット
fprintf('\n実験完了。プロットを作成中...\n');

total_samples = num_sets * num_ir * 3; 
plot_SDR = zeros(total_samples, 1);
plot_IR  = cell(total_samples, 1);
plot_Method = cell(total_samples, 1);

unique_labels = cell(1, num_ir);
for ir_idx = 1:num_ir
    unique_labels{ir_idx} = sprintf('IR %d (%s)', ir_idx, rt60_labels{ir_idx});
end

idx = 1;
for ir_idx = 1:num_ir
    for s_idx = 1:num_sets
        % CR法
        plot_SDR(idx) = delta_SDR_CR(s_idx, ir_idx);
        plot_IR{idx}  = unique_labels{ir_idx};
        plot_Method{idx} = 'Proposed CR';
        idx = idx + 1;
        
        % Lifting法
        plot_SDR(idx) = delta_SDR_Lift(s_idx, ir_idx);
        plot_IR{idx}  = unique_labels{ir_idx};
        plot_Method{idx} = 'Lifting';
        idx = idx + 1;
        
        % WPE法
        plot_SDR(idx) = delta_SDR_WPE(s_idx, ir_idx);
        plot_IR{idx}  = unique_labels{ir_idx};
        plot_Method{idx} = 'Monaural WPE';
        idx = idx + 1;
    end
end

plot_IR = categorical(plot_IR, unique_labels);
plot_Method = categorical(plot_Method, {'Proposed CR', 'Lifting', 'Monaural WPE'});

figure('Position', [100, 100, 750, 480]);
boxchart(plot_IR, plot_SDR, 'GroupByColor', plot_Method);

grid on;
ylabel('\Delta SDR [dB]', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Room Impulse Response (RIR)', 'FontSize', 12);
title('Dereverberation Performance Across Different RT60', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 10);
set(gca, 'FontSize', 11);

fprintf('\n---- 平均処理時間 (1セットあたり) ----\n')
fprintf('Proposed CR : %.4f 秒\n', mean(time_CR(:)));
fprintf('Lifting     : %.4f 秒\n', mean(time_Lift(:)));
fprintf('Monaural WPE: %.4f 秒\n', mean(time_WPE(:)));