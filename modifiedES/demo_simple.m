s1 = audioread("../data/ViolaL20210803_8000Hz_1.wav");
%s1 = s1(1:fs); soundsc(s1, fs); pause(2)
[s2, fs] = audioread("../data/ViolaL20210803_8000Hz_2.wav");
%s2 = s2(1:fs); soundsc(s2, fs);

[ir, fsir] = audioread("../data/BF TL SPACE LIBRARY/Drumbrella/Drumbrella 5'.L.wav");
ir = ir(1:reverb_time(ir,fsir)*fsir);
% [ir, fsir] = audioread("../data/BF TL SPACE LIBRARY/Drumbrella/Drumbrella 7'.R.wav");
% ir = ir(1:8000);
ir8000 = resample(ir, fs,fsir);
lh = length(ir8000);

x1 = conv(s1, ir8000);
x2 = conv(s2, ir8000);
l1 = length(s1); l2 = length(s2);
soundsc(s1, fs); pause(length(s1)/fs+0.5); soundsc(x1, fs)
%% blind deconvolution

S = [convmtx(x2, l1), convmtx(-x1, l2)];
[~,~,V] = svd(S, "econ");
shat = V(:,end);
shat1 = shat(1:l1);
shat2 = shat(l1+1:end);
A = [convmtx(shat1, lh); convmtx(shat2, lh)];
b = [x1; x2];
hhat = A \ b;

%% 
long_s = audioread("../data/brahms_vc_8000Hz.wav");
xx = conv(long_s, ir8000);
XX = fft(xx);
hhat_ext = zeros(size(xx)); 
hhat_ext(1:lh) = hhat;
Hhat = fft(hhat_ext);
z = ifft(XX ./ Hhat);
soundsc(xx, fs); 
pause(round(length(xx)/fs));
soundsc(z, fs)

% Spectral Distance
LogSpectralDistance(long_s, xx, fs)
LogSpectralDistance(long_s, z, fs)

%% 

tiledlayout(2,1)
nexttile;
lin = linspace(0,2*pi);
r = roots(hhat);
scatter(real(r), imag(r), "DisplayName", "roots of hhat")
hold on;
plot(cos(lin), sin(lin), "DisplayName", "unitcircle")
xlabel("real part"); ylabel("imag part")
legend

nexttile;
rr = roots(ir8000);
scatter(real(rr), imag(rr), "DisplayName", "roots of ir");
hold on
plot(cos(lin), sin(lin), "DisplayName", "unitcircle")
xlabel("real part"); ylabel("imag part"); legend

% どちらも単位円の内外に〇点を持っている