[long_s, fs] = audioread("../data/brahms_vc_8000Hz.wav");
[ir, fsir] = audioread("../data/BF TL SPACE LIBRARY/Drumbrella/Drumbrella 5'.L.wav");
ir = ir(1:ceil(reverb_time(ir,fsir)*fsir));
ir8000 = resample(ir, fs, fsir);
lh = length(ir8000);

y = conv(long_s, ir8000);
h = ir8000;

%% convenional deconv
% recovering 1
h_ext = zeros(size(y)); 
h_ext(1:lh) = h;
H1 = fft(h_ext);
z1 = ifft( fft(y) ./ H1);

% recovering 2
M = convmtx(h, length(h));
e = zeros(size(M,1), 1); e(1) = 1;
h_inv = M \ e;
z2 = conv(y, h_inv);

%% oppenheim
P = length(h); 
L = 1024*8;  % frame length
Q = 0;     % hamidasi
Lp = P -1 + Q;
La = L + Lp;

% H_L = fft(h, L);
% % 正則化と閾値処理
% if options.adaptiveThreshold
%     epsilon = 0.01 * max(abs(H_L));  % 動的設定
% else
%     epsilon = options.epsilon;       % 固定値
% end
% 
% if options.regularization
%     % 高度な正則化（Tikhonov型）
%     H_L_safe = conj(H_L) ./ (abs(H_L).^2 + epsilon);
% else
%     % 単純なしきい値処理
%     H_L_safe = H_L;
%     H_L_safe(abs(H_L) < epsilon) = epsilon;
% end
% R = ceil(length(y)/(La));
% d  = zeros(L * R, 1);
% y_ext = zeros(R * La, 1);
% y_ext(1:length(y)) = y;
% 
% n_range = 1:L;
% for r = 0:R-1
%     % #r frame proc.
%     if r < R-1
%         yhat_r = y_ext(n_range + r * La) + y_ext(n_range + r * La + L);
%     else
%         yhat_r = y_ext(n_range + r * La);
%     end
%     y_rL = yhat_r; y_rL(L+1:end) = 0;
%     Y_rL = fft(y_rL, L);
%     if options.regularization
%         D_r = Y_rL .* H_L_safe;
%     else
%         D_r = Y_rL ./ H_L_safe;
%     end
%     dr = ifft(D_r, L);
%     d(n_range + r*La ) = dr;
% end


d = efficient_deconvolution(y,h,L,Q,options);

%% analyzing the results
% listening
% soundsc(y, fs); 
% pause(round(length(y)/fs));
% soundsc(z1, fs)
% pause(round(length(y)/fs));
% soundsc(z2, fs)
% pause(round(length(y)/fs));
% soundsc(d, fs)

% Spectral Distance
LogSpectralDistance(long_s, y, fs)
LogSpectralDistance(long_s, z1, fs)
LogSpectralDistance(long_s, z2, fs)
LogSpectralDistance(long_s, d, fs)


%
figure; tiledlayout(5,1);
nexttile; plot(long_s); xlabel("dry"); xlim([1,length(y)])
nexttile; plot(y); xlabel("wet"); xlim([1,length(y)])
nexttile; plot(z1); xlabel("inverse in the FFT domain"); xlim([1,length(y)])
nexttile; plot(z2); xlabel("inverse using filter"); xlim([1,length(y)])
nexttile; plot(d); xlabel("oppenheim"); xlim([1,length(y)])
