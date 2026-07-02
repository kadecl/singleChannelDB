% 1. ダミーデータ作成（手元にある青い波形データの代わり）
% ※実際にお手元のデータを読み込む際は、以下の3行を[y, fs] = audioread('filename.wav'); 等に置き換えてください。
[y, fs] = audioread("/Users/itsukiyashima/Documents/MATLAB/singleChannelDB/data/brahms_vc.wav");
y = y(70000:end);
y(90001:130000) = y(90001:130000) .* (0.9999).^(1:40000)';
% --- ここから黒いアイコン風波形への変換処理 ---

% 2. 上下の最大値を結ぶ包絡線（エンベロープ）の計算
[y_upper, ~] = envelope(y, 300, 'peak'); % 窓幅300で滑らかに抽出

% 3. 概念図として見栄えが良い点数（数〜数百点）にダウンサンプリング
% アイコン風の「縦棒」の密度を調整するパラメータです
num_bars = 30; 
ds_factor = floor(length(y_upper) / num_bars);
y_ds = y_upper(1:ds_factor:end);

% 4. 描画用の擬似的な縦棒データの作成 (stem関数を応用)
x_bars = 1:length(y_ds);

% 5. 描画設定（背景透明・枠線なしのプロ仕様）
fig = figure('Color', 'none'); % Figureの背景を透明に設定
% fig = figure;
ax = axes('Color', 'none');    % Axes（軸）の背景も透明に設定

% 縦棒（黒）として描画
h = stem(x_bars, y_ds, 'LineWidth', 10, 'Color', 'k', 'Marker', 'none');
hold on;
% 反転した下側の縦棒も描画して上下対称の波形アイコンにする
stem(x_bars, -y_ds, 'LineWidth', 10, 'Color', 'k', 'Marker', 'none');

% 軸の目盛りや枠線を完全に消去
axis off; 

% 6. 背景を透明な状態のままPNGとして保存する小技
set(fig, 'InvertHardcopy', 'off'); % 保存時に背景を白に強制反転させない設定
exportgraphics(fig, 'waveform_icon.png', 'BackgroundColor', 'none', 'Resolution', 300);

disp('背景透過の黒い波形画像（waveform_icon.png）を出力しました！');