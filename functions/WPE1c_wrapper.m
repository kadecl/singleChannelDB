function y_out = WPE1c_wrapper(x_in, fs, F, Lir)
% WPE_MONAURAL_WRAPPER  NTT製PコードWPEを単一チャネル(モノラル)信号に対して実行する
%
% 入力:
%   x_in : モノラル観測信号ベクトル [サンプル数 x 1] 
%   fs   : サンプリング周波数 (Hz)
% 出力:
%   y_out: 残響除去後の信号ベクトル [サンプル数 x 1]

    % 1. 入力ベクトルの形状チェック (必ず [サンプル数 x 1] の縦ベクトルにする)
    if size(x_in, 2) > 1
        x_in = x_in.';
    end
    
    % 2. 設定ファイル名を確定 (型エラーを防ぐため完全なchar型にする)
    % wpe.pが認識しやすいよう、余計なサブフォルダを作らずカレントディレクトリに配置します
    cfg_filepath = char(fullfile(pwd, 'wpe_local_monaural_temp.m'));
    
    % 3. モノラル用にパラメータをチューニングした設定ファイルを動的に生成
    fid = fopen(cfg_filepath, 'w');
    if fid == -1
        error('一時設定ファイルの作成に失敗しました。');
    end
    
    fprintf(fid, '%% モノラル実験用に動的生成された設定ファイル\n');
    fprintf(fid, 'fs = %d;\n', fs);
    fprintf(fid, 'num_mic = 1;\n');        % モノラル入力に固定
    fprintf(fid, 'num_out = 1;\n');        % モノラル出力に固定
    fprintf(fid, 'blk_len = 100;\n');      % 発話バッチ処理用
    fprintf(fid, 'opt_blk_sz = 1;\n');
    
    % STFTパラメータ
    fprintf(fid, 'analy_param = struct(''win_size'', %d, ''shift_size'', %d, ''win'', hanning(512));\n', F.FFTnum, F.shift);
    
    % 予測フィルタの設定 [タップ数; 予測遅延; 上限周波数]
    fprintf(fid, 'channel_setup = [%d; 3; inf];\n', Lir); 
    
    % 推定コンフィグ (p_channelはモノラルなので 1 のみ)
    fprintf(fid, 'ssd_param = struct(''channel_setup'', channel_setup, ''p_channel'', 1, ''speech_order'', 20);\n');
    fprintf(fid, 'ssd_conf = struct(''max_iter'', 3, ''spcorr'', ''scaleye'', ''scaling'', 1, ''forget'', 0.7);\n');
    fprintf(fid, 'enh_conf = struct(''method'', ''lti'', ''osub'', 1.0, ''scal'', 2, ''floor'', -80);\n');
    
    fclose(fid);
    
    % 4. 暗号化されているNTT公式WPE(Pコード)を実行
    % 引数のパスが確実なテキストスカラー(char)であることを保証する
    try
        y_out = wpe(x_in, cfg_filepath);
    catch ME
        if exist(cfg_filepath, 'file'), delete(cfg_filepath); end
        rethrow(ME);
    end
    
    % 5. 使用した一時ファイルのクリーンアップ
    if exist(cfg_filepath, 'file')
        delete(cfg_filepath);
    end
end