%%% 20250228 最も基礎的なあいでぃあ．長めのリバーブに対してシングルチャネルでできているので新規性はアピールできそう．
%%% 一方で，


clear all
%% settings
% DGT tool
F = DGTtool("windowLength", 512, "windowShift", 128);

% audio data loading
num_piece = 2;
for i=1:num_piece
    audiofilename = "input/9w3r7y_01_%d.wav";
    [temp, fs] =  audioread(sprintf(audiofilename, i));
    ss{i} = temp(:,1);
end
s = [ss{1}; ss{2}];
num_mic = 1;

% ir data loading
[ir, ~] = audioread("BFTLSPACELIBRARY_ChamberAL_11025.wav");
%[ir, fs] = audioread("logic03sWoodenBooth8000.wav");
%[ir{2}, fs(2)] = audioread("../data/BF TL SPACE LIBRARY/Drumbrella/Drumbrella 5'.R.wav");

%% slra options
useSLRA = 0;
useAmplitude = 1;
avoidNearZeroPolys = 1;
opt.MaxIterations = 100;
opt.epsgrad = 1e-12;
opt.solver = "c";

%% malloc
rt = zeros(num_mic, 1); % reverb time in sec
obs = cell(num_mic,num_piece); % obs
X = cell(num_mic, num_piece); % STFTs of obs. signal
N = zeros(num_mic, num_piece); % length of obserbation
L = zeros(num_mic, 1); % ir length = d + 1
l = zeros(num_mic, num_piece); % length of deconvoluted signal

%% resample and shortening
for i = 1:num_mic
    % calcurate the reververation time RT60
    rt(i) = reverb_time(ir, fs);
    ir = ir(1:fs); ir = ir / norm(ir);
    L(i) = size(F(ir), 2);

    for j = 1:num_piece
        obs{i, j} = conv(ss{j}, ir);
        X{i, j} = F(obs{i, j});
        N(i, j) = size(X{i, j}, 2);
    end
end
obsCat = [obs{1}; obs{2}];
d = L -1;

%% deconvolution in TF domain
% malloc
Fq = size(X{1}, 1);
for i = 1:2, l(i) = N(1, i) - L(1) + 1; end
RET1 = zeros(Fq, l(1)); RET2 = zeros(Fq, l(2));
RET_slra = zeros(Fq, sum(l));
% minimum to 2nd minimum singular value ratio
singValRatio = zeros(Fq, 1);

% frequency-wise processing
parfor f = 1:Fq
    temp = X{1, 1};
    X1f = temp(f, :)';
    temp = X{1, 2};
    X2f = temp(f, :)';
    Xf = [X1f; X2f];
    XfNorm = norm(Xf);
    XfPhase = Xf ./ abs(Xf);

    % (not) structured low-rank approximation
    Syl = [convmtx(X2f, l(1)) convmtx(-X1f, l(2))];
    [~, Sigma, V] = svd(Syl, "econ");
    Sigma = diag(Sigma); singValRatio(f) = Sigma(end) / Sigma(end-1);
    Vf = V(:,end);
    Vf = XfNorm * Vf / norm(Vf); % scale recovery
    v1 = Vf(1:l(1)); 
    v2 = Vf(l(1)+1:end);

    % phase normalization
    A = [convmtx(v1, L(1)); convmtx(v2, L(1))];
    Hhat = A \ Xf; % calc. IR (TF domain)
    HhatPhase = Hhat ./ abs(Hhat);
    [Hhatmax, idx_HhatPeak] = max(abs(Hhat));
    v1 = v1 * HhatPhase(idx_HhatPeak); v2 = v2 * HhatPhase(idx_HhatPeak); 
    RET1(f,:) = v1';   RET2(f,:) = v2';

    % if useSLRA % && avoidNearZeroPolys
    %     if useAmplitude
    %         A1f = abs(X1f); A2f = abs(X2f);
    %         [ph, info] = gcd_nls({A1f, A2f}, [], d, opt);
    %         Hhat = info.Rh(:);
    %     else
    %         % どうやら，カーネルが二次元なのでかなり時間がかかる
    %         [ph, info] = gcd_nls_complex({X1f, X2f}, [], d, opt);
    %         Hhat = info.hh(:);
    %     end
    %     Syl = fullSyl(ph,d);
    %     [U, Sigma, V] = svd(Syl, "econ");
    %     Vtemp = V(:,end);
    %     v1 = Vtemp(1:l(1)); v1 = a1 * v1 / norm(v1); % scale recovery
    %     v2 = Vtemp(l(1)+1:end); v2 = a2 * v2 / norm(v2); % scale recovery
    % 
    %     [HhatPeak, idx_HhatPeak] = max(abs(Hhat));
    %     v1 = v1 / Hhat(idx_HhatPeak); v2 = v2 / Hhat(idx_HhatPeak); % phase normalization
    %     RET_slra(1,f,:) = v1';  RET_slra(2,f,:) = v2';
    % 
    %     if useAmplitude
    %         % 単純なLRAで計算した際の位相を貼り付ける
    %         RET_f = RET(f,:);
    %         RET_f_phase = RET_f / abs(RET_f);
    %         RET_slra(f,:) = RET_slra(f,:) .* RET_f_phase;
    %     end
    % end
end
%% 
ret = cell(1);
ret{1} = F.pinv(RET1(:,:));  ret{2} = F.pinv(RET2(:,:));
for i=1:2
    
    % ret_slra = F.pinv(RET_slra);
    len_ss = length(ss{i});
    [SDR,SIR,SAR,perm] = bss_eval_sources(obs{i}(1:len_ss)', ss{i}');
    fprintf("obs: SDR %2.3f, SIR %2.3f, SAR %2.3f\n", SDR, SIR, SAR)

    [SDR,SIR,SAR,perm] = bss_eval_sources(ret{i}(1:len_ss)', ss{i}');
    fprintf("ret: SDR %2.3f, SIR %2.3f, SAR %2.3f\n", SDR, SIR, SAR)

    audiowrite(sprintf("output/wet_%d.wav", i), obs{i}, fs)
    audiowrite(sprintf("output/lra_%d.wav", i), ret{i}, fs)
    % audiowrite("output/slra.wav", ret_slra, fs)
end

F.plotReassign(ss{1});title("dry"); saveas(gcf, audiofilename + "dry.png")
F.plotReassign(obs{1}); title("wet"); saveas(gcf, audiofilename + "wet.png")
F.plotReassign(ret{1}); title("lra"); saveas(gcf, audiofilename + "ret.png")
if useSLRA, F.plotReassign(ret_slra); title("slra");end

% obs: SDR 0.689, SIR Inf, SAR 0.689
% ret: SDR 4.171, SIR Inf, SAR 4.171
% obs: SDR -0.120, SIR Inf, SAR -0.120
% ret: SDR 4.543, SIR Inf, SAR 4.543
% audiowrite("output/slra.wav", ret_slra, fs)
% player = audioplayer(obsCat, 8000, 16, 4); play(player);

