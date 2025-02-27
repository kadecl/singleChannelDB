[ir, fsir] = audioread("../data/BF TL SPACE LIBRARY/Drumbrella/Drumbrella 5'.L.wav");
ir = ir(1:ceil(reverb_time(ir,fsir)*fsir));
ir8000 = resample(ir, fs,fsir);
lh = length(ir8000);

long_s = audioread("../data/brahms_vc_8000Hz.wav");
y = conv(long_s, ir8000);

h = ir8000;
P = length(h);
L = 2*P;
Lp = P-1;
La = L + Lp;

H = fft(h, L);
R = ceil(length(y)/(La));
d  = zeros(L * R, 1);
y_ext = zeros(R * La, 1);
y_ext(1:length(y)) = y;

n_range = 1:L;
for r = 0:R-1
    
    yrL = y_ext(n_range + r * La);
    YrL = fft(yrL, L);
    dr = ifft(YrL ./ H, L);
    d(n_range + r*L ) = dr;
end

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
LogSpectralDistance(long_s, z, fs)
LogSpectralDistance(long_s, z2, fs)
LogSpectralDistance(long_s, d, fs)


%
figure; tiledlayout(5,1);
nexttile; plot(long_s); xlabel("dry"); xlim([1,length(y)])
nexttile; plot(y); xlabel("wet"); xlim([1,length(y)])
nexttile; plot(z1); xlabel("inverse in the FFT domain"); xlim([1,length(y)])
nexttile; plot(z2); xlabel("inverse using filter"); xlim([1,length(y)])
nexttile; plot(d); xlabel("oppenheim"); xlim([1,length(y)])
