%%% 20250228 demo.m に，SylvesterSLRAを使う.
%%% 位相のオラクルを使う

clear all
%% settings
% DGT tool
F = DGTtool("windowLength", 512, "windowShift", 128);

% audio data loading
num_piece = 2;
for i=1:num_piece
    audiofilename = "9w3r7y_01_%d";
    [temp, fs] =  audioread("input/" + sprintf(audiofilename,i)+".wav");
    ss{i} = temp(:,1);
end
s = [ss{1}; ss{2}];
num_mic = 1;

% ir data loading
[ir, ~] = audioread("BFTLSPACELIBRARY_ChamberAL_11025.wav");
%[ir, fs] = audioread("logic03sWoodenBooth8000.wav");
%[ir{2}, fs(2)] = audioread("../data/BF TL SPACE LIBRARY/Drumbrella/Drumbrella 5'.R.wav");

%% slra options
useSLRA = 1;
useAmplitude = 1;
avoidNearZeroPolys = 1;
opt.MaxIterations = 100;
opt.epsgrad = 1e-12;
opt.solver = "c";

%% malloc
rt = zeros(num_mic, 1); % reverb time in sec
obs = cell(num_mic,num_piece); % obs
X = cell(num_mic, num_piece); % STFTs of obs. signal
S = cell(num_mic, num_piece); % STFTs of dry signal
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
        S{i, j} = F(ss{j});
    end
end
obsCat = [obs{1}; obs{2}];
d = L -1;

%% deconvolution in TF domain
% malloc
Fq = size(X{1}, 1);
for i = 1:2, l(i) = N(1, i) - L(1) + 1; end
RET1 = zeros(Fq, l(1)); RET2 = zeros(Fq, l(2));
RET_slra1 = zeros(Fq, l(1));  RET_slra2 = zeros(Fq, l(2));
% minimum to 2nd minimum singular value ratio
singValRatio = zeros(Fq, 1);

% frequency-wise processing
parfor f = 1:Fq
    temp = X{1, 1};
    X1f = temp(f, :)';
    % a1 = norm(X1f); X1fとX2fの大きさで別々に正規化してたらだめじゃん．
    temp = X{1, 2};
    X2f = temp(f, :)';
    % a2 = norm(X2f);
    Xf = [X1f; X2f];
    XfNorm = norm(Xf);
    XfPhase = Xf ./ abs(Xf);

    % (not) structured low-rank approximation
    Syl = [convmtx(X2f, l(1)) convmtx(-X1f, l(2))];
    [Uf, Sigmaf, Vf] = svd(Syl, "econ");
    Sigma = diag(Sigmaf); singValRatio(f) = Sigma(end) / Sigma(end-1);
    Vfend = Vf(:,end);
    Vfend = XfNorm * Vfend / norm(Vfend); % scale recovery
    VfAng = angle(Vfend);

    % % IR estimation
    % A = [convmtx(v1, L(1)); convmtx(v2, L(1))];     b = [X1f; X2f];
    % Hhat = A \ b; % calc. IR (TF domain)
    % HhatPhase = Hhat / abs(Hhat);
    % [Hhatmax, idx_HhatPeak] = max(abs(Hhat));

    % use Oracle PHASE
    S1f = S{1,1}(f,:); S2f = S{1,2}(f,:);  Sf = [S1f'; S2f'];
    SfAng = angle(Sf);
    % aveAngDiff = mean(SfAng - VfAng);
    weight = abs(Vfend);
    aveAngDiff = sum(weight.*(SfAng - VfAng)) / sum(weight);
    EstPhaseUsingOracle = exp(1i * aveAngDiff);
    Vfend = Vfend .* EstPhaseUsingOracle; % / Hhatmax;

    % v1 = v1 / Hhatmax; v2 = v2 / Hhatmax; ::
    % ただでさえあやしい時間周波数領域での畳み込みの成立を前提とした
    % 誤差が蓄積しているHhatを使って正規化するのはどうなん．．ということでつかわないかも．
    v1 = Vfend(1:l(1));
    v2 = Vfend(l(1)+1:end);
    RET1(f,:) = v1';   RET2(f,:) = v2';

    if useSLRA
        % [ph, info] = gcd_syl({X1f, X2f}, [], d, opt); owannneeeee
        % Hhat = info.hh(:);
        ph = {X1f, X2f};
        for iter = 1:5
            Sigmaf(end, end) = 0;
            ph = orthoProjOntoSylvester(Uf * Sigmaf * Vf', N-1, d);
            [Uf, Sigmaf, Vf] = svd(fullSyl(ph,d), "econ");
        end
        Vfend = Vf(:,end);
        Vfend = XfNorm * Vfend / norm(Vfend); % scale recovery
        VfAng = angle(Vfend);

        aveAngDiff = mean(SfAng - VfAng);
        EstPhaseUsingOracle = exp(1i * aveAngDiff);
        Vfend = Vfend .* EstPhaseUsingOracle;

        RET_slra1(f,:) = Vfend(1:l(1))';  RET_slra2(f,:) = Vfend(l(1)+1:end)';
    end
end
%% 
ret = cell(1);
ret{1} = F.pinv(RET1(:,:));  ret{2} = F.pinv(RET2(:,:));
ret_slra{1} = F.pinv(RET_slra1(:,:));  ret_slra{2} = F.pinv(RET_slra2(:,:));
for i=1:2
    
    % ret_slra = F.pinv(RET_slra);
    len_ss = length(ss{i});
    [SDR,SIR,SAR,perm] = bss_eval_sources(obs{i}(1:len_ss)', ss{i}');
    fprintf("obs: SDR %2.3f, SIR %2.3f, SAR %2.3f\n", SDR, SIR, SAR)

    [SDR,SIR,SAR,perm] = bss_eval_sources(ret{i}(1:len_ss)', ss{i}');
    fprintf("ret: SDR %2.3f, SIR %2.3f, SAR %2.3f\n", SDR, SIR, SAR)

    [SDR,SIR,SAR,perm] = bss_eval_sources(ret_slra{i}(1:len_ss)', ss{i}');
    fprintf("SLRA: SDR %2.3f, SIR %2.3f, SAR %2.3f\n", SDR, SIR, SAR)

    audiowrite("output/" + sprintf(audiofilename +"_wet_ORACLE.wav", i), obs{i}, fs)
    audiowrite("output/" + sprintf(audiofilename +"_lra_ORACLE.wav", i), ret{i}, fs)
    audiowrite("output/" + sprintf(audiofilename +"_slra_ORACLE.wav", i), ret_slra{i}, fs)
end

F.plotReassign(ss{1});title("dry"); 
saveas(gcf, "result/" + sprintf(audiofilename,1) + "_dry.png")
F.plotReassign(obs{1}); title("wet"); 
saveas(gcf, "result/" + sprintf(audiofilename,1) + "_wet.png")
F.plotReassign(ret{1}); title("lra");
saveas(gcf, "result/" + sprintf(audiofilename,1) + "_LRA_ORACLE.png")
if useSLRA
    F.plotReassign(ret_slra{1}); title("slra");
    saveas(gcf, "result/" + sprintf(audiofilename,1) + "_SLRA_ORACLE.png")
end

%どの結果も時間のマイナス方向に信号が漏れ出ているのがきになる

% audiowrite("output/slra.wav", ret_slra, fs)
% player = audioplayer(obsCat, 8000, 16, 4); play(player);

% obs: SDR 0.689, SIR Inf, SAR 0.689
% ret: SDR 6.214, SIR Inf, SAR 6.214
% SLRA: SDR -14.830, SIR Inf, SAR -14.830
% obs: SDR -0.120, SIR Inf, SAR -0.120
% ret: SDR 5.250, SIR Inf, SAR 5.250
% SLRA: SDR -15.199, SIR Inf, SAR -15.199