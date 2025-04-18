%%% 20250228 demo.m に，オラクル情報を使った差分
%%% 残響除去後のスペクトルに，オラクルから得る位相差分を貼り付ける方法を検証する．位相の話がいったん落ち着くので，
%%% 手法のメインなアイディアの良さを評価できることが期待できる．
%%% 松山さんの直行射影のアイディアをとりいれる。
%%% TF領域でのサブバンドごとのAGCD計算により得たインパルス応答を直行射影してIRを復元する。そのあと再度クリーンを計算する

clear
%% loading
% audio data loading
num_piece = 2;
for i=1:num_piece
    audiofilename = "9w3r7y_01_%d";
    [temp, fs] =  audioread("input/" + sprintf(audiofilename,i)+".wav");
    s{i} = temp(:,1);
end
num_mic = 1;

% ir data loading
[ir, ~] = audioread("BFTLSPACELIBRARY_ChamberAL_11025.wav");
%[ir, fs] = audioread("logic03sWoodenBooth8000.wav");
%[ir{2}, fs(2)] = audioread("../data/BF TL SPACE LIBRARY/Drumbrella/Drumbrella 5'.R.wav");

%% malloc
obs = cell(num_mic, num_piece); % obs
X = cell(num_mic, num_piece); % STFTs of obs. signal
S = cell(num_mic, num_piece);   % STFTs of dry signal
N = zeros(num_mic, num_piece); % length of obserbation
Lir = zeros(num_mic, 1); % ir length = d + 1
l = zeros(num_mic, num_piece); % length of deconvoluted signal

%% convolution
% calcurate the reververation time RT60 (IR shortening)
rt = reverb_time(ir, fs);
ir = ir(1:fs); ir = ir / norm(ir);
IR = dgt(ir,gt,a,M,'timeinv');
Lir = size(IR, 2);

for j=1:num_piece
    obs{j} = postpad(obs{j}, n);
    X{j} = dgt(obs{j},gt,a,M,'timeinv');
    N(j) = size(X{j}, 2);
end
obsCat = [obs{1}; obs{2}];

%% resample and shortening
d = Lir - 1;
padding = 0;
d = d - padding;
Lir = Lir - padding;

%% deconvolution in TF domain
Fq = size(X{1}, 1);
l = N - Lir + 1;

RET1 = zeros(Fq, l(1)); RET2 = zeros(Fq, l(2));
RET1W = zeros(Fq, l(1)); RET2W = zeros(Fq, l(2));
% minimum to 2nd minimum singular value ratio
singValRatio = zeros(Fq, 1);

% frequency-wise processing
for f = 1:Fq
    temp = X{1, 1};  X1f = temp(f, :)';
    temp = X{1, 2};   X2f = temp(f, :)';
    Xf = [X1f; X2f];
    XfNorm = norm(Xf);
    XfPhase = Xf ./ abs(Xf);

    % (not) structured low-rank approximation
    Syl = [convmtx(X2f, l(1)) convmtx(-X1f, l(2))];
    [~, Sigma, V] = svd(Syl, "econ");
    Sigma = diag(Sigma); singValRatio(f) = Sigma(end) / Sigma(end-1);
    Vf = V(:,end);
    Vf = XfNorm * Vf / norm(Vf); % scale recovery
    VfAng = angle(Vf);

    % % IR estimation
    % A = [convmtx(v1, L(1)); convmtx(v2, L(1))];     b = [X1f; X2f];
    % Hhat = A \ b; % calc. IR (TF domain)
    % HhatPhase = Hhat / abs(Hhat);
    % [Hhatmax, idx_HhatPeak] = max(abs(Hhat));

    % use Oracle PHASE
    S1f = S{1,1}(f,:); S2f = S{1,2}(f,:);  Sf = [S1f'; S2f'];
    SfAng = angle(Sf);
    %aveAngDiff = mean(SfAng - VfAng);
    aveAngDiff = weightedMeanWithZeroPadding(SfAng, VfAng);
    EstPhaseUsingOracle = exp(1i * aveAngDiff);
    Vf_avephase = Vf * EstPhaseUsingOracle; % / Hhatmax;

    weight = abs(Sf);
    % aveAngDiffW = sum(weight.*(SfAng - VfAng)) / sum(weight);
    aveAngDiffW = weightedMeanWithZeroPadding(SfAng, VfAng, weight);
    EstPhaseUsingOracleW = exp(1i * aveAngDiffW);
    Vf_w = Vf * EstPhaseUsingOracleW;

    % v1 = v1 / Hhatmax; v2 = v2 / Hhatmax; ::
    % ただでさえあやしい時間周波数領域での畳み込みの成立を前提とした
    % 誤差が蓄積しているHhatを使って正規化するのはどうなん．．ということでつかわないかも．
    v1 = Vf_avephase(1:l(1));% v1 = a1 * v1 / norm(v1); % scale recovery
    v2 = Vf_avephase(l(1)+1:end); %v2 = a2 * v2 / norm(v2); % scale recovery
    RET1(f,:) = v1';   RET2(f,:) = v2';

    RET1W(f,:) = Vf_w(1:l(1))'; RET2W(f,:) = Vf_w(l(1)+1:end)';

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
ret_w{1} = F.pinv(RET1W); ret_w{2} = F.pinv(RET2W);

%%%===
% orthogonal projection 
ltfatstart
gt=gabwin('hann',a,M);
corwin=xcorr(circshift(gt,M/2),circshift(gd,M/2));
K=size(Gfull,2)*a; % K : IR length
%corwin= buffer(corwin,K);
[~,I]=(max(corwin,[],"all"));
corwin=circshift(corwin,-I+1);  % window for B
cd = gabdual(corwin,a,M);   % window for B_inv

hre=idgtreal(G,cd,a,M,'timeinv');  % G:射影前のフィルタ
J = dgt(hre,gd,1,M,'timeinv');
Hobl = dgt(permute(J,[2 1]),gt,a,M,'timeinv');
%%%===

for i=1:2
    
    % ret_slra = F.pinv(RET_slra);
    len_ss = length(s{i});
    len_ret = length(ret{i});
    if len_ss > len_ret
        s{i} = s{i}(1:len_ret);
        len_ss = len_ret;
    end
    [SDR,SIR,SAR,perm] = bss_eval_sources(obs{i}(1:len_ss)', s{i}');
    fprintf("obs: SDR %2.3f, SIR %2.3f, SAR %2.3f\n", SDR, SIR, SAR)

    [SDR,SIR,SAR,perm] = bss_eval_sources(ret{i}(1:len_ss)', s{i}');
    fprintf("ret: SDR %2.3f, SIR %2.3f, SAR %2.3f\n", SDR, SIR, SAR)

    [SDR,SIR,SAR,perm] = bss_eval_sources(ret_w{i}(1:len_ss)', s{i}');
    fprintf("ret(w): SDR %2.3f, SIR %2.3f, SAR %2.3f\n", SDR, SIR, SAR)

    audiowrite("output/" + sprintf(audiofilename +"_wet.wav", i), obs{i}, fs)
    audiowrite("output/" + sprintf(audiofilename +"_lra_ORACLE_pad=%d.wav", i, padding), ret{i}, fs)
    audiowrite("output/" + sprintf(audiofilename +"_lra_ORACLE(w)_pad=%d.wav", i, padding), ret_w{i}, fs)
    % audiowrite("output/slra.wav", ret_slra, fs)
end

F.plotReassign(s{1});title("dry"); 
saveas(gcf, "result/" + sprintf(audiofilename,1) + "_dry.png")
F.plotReassign(obs{1}); title("wet"); 
saveas(gcf, "result/" + sprintf(audiofilename,1) + "_wet.png")
F.plotReassign(ret{1}); title("lra");
saveas(gcf, "result/" + sprintf(audiofilename + "_ORACLE_ret_pad=%d.png",1, padding) )
F.plotReassign(ret_w{1}); title("lra(w)");
saveas(gcf, "result/" + sprintf(audiofilename + "_ORACLE(w)_ret.png",1, padding) )
if useSLRA, F.plotReassign(ret_slra); title("slra");end


% audiowrite("output/slra.wav", ret_slra, fs)
% player = audioplayer(obsCat, 8000, 16, 4); play(player);
%% non-weighted mean  
% obs: SDR 0.689, SIR Inf, SAR 0.689
% ret: SDR 6.136, SIR Inf, SAR 6.136
% ret(w): SDR 6.214, SIR Inf, SAR 6.214
% obs: SDR -0.120, SIR Inf, SAR -0.120
% ret: SDR 5.033, SIR Inf, SAR 5.033
% ret(w): SDR 5.250, SIR Inf, SAR 5.250
% さらに，位相差分の重みをオラクルの振幅を使用するとわずかに改善
% ret(w): SDR 6.296, SIR Inf, SAR 6.296
% ret(w): SDR 5.102, SIR Inf, SAR 5.102


%どの結果も時間のマイナス方向に信号が漏れ出ているのがきになる
% weight をいれることで若干改善．ただ，きいてもわからない

% padding =
% 
%     10
% 
% index computed (OLA)
% obs: SDR 0.689, SIR Inf, SAR 0.689
% ret: SDR -0.123, SIR Inf, SAR -0.123
% ret(w): SDR -0.151, SIR Inf, SAR -0.151
% obs: SDR -0.120, SIR Inf, SAR -0.120
% ret: SDR 3.056, SIR Inf, SAR 3.056
% ret(w): SDR 2.829, SIR Inf, SAR 2.829
% index computed (OLA)
% 2
% obs: SDR 0.689, SIR Inf, SAR 0.689
% ret: SDR 8.056, SIR Inf, SAR 8.056
% ret(w): SDR 8.071, SIR Inf, SAR 8.071
% obs: SDR -0.120, SIR Inf, SAR -0.120
% ret: SDR 7.384, SIR Inf, SAR 7.384
% ret(w): SDR 7.209, SIR Inf, SAR 7.209
% -2 (じっさいよりながいIR長)
% index computed (OLA)
% obs: SDR 0.668, SIR Inf, SAR 0.668
% ret: SDR 4.215, SIR Inf, SAR 4.215
% ret(w): SDR 4.288, SIR Inf, SAR 4.288
% obs: SDR -0.136, SIR Inf, SAR -0.136
% ret: SDR 2.163, SIR Inf, SAR 2.163
% ret(w): SDR 2.054, SIR Inf, SAR 2.054
% -10
% index computed (OLA)
% obs: SDR 0.825, SIR Inf, SAR 0.825
% ret: SDR -1.172, SIR Inf, SAR -1.172
% ret(w): SDR -1.099, SIR Inf, SAR -1.099
% obs: SDR -0.028, SIR Inf, SAR -0.028
% ret: SDR -1.837, SIR Inf, SAR -1.837
% ret(w): SDR -1.775, SIR Inf, SAR -1.775