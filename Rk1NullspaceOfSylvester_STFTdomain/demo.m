clear all
%% settings
% DGT tool
F = DGTtool("windowLength", 512, "windowShift", 128);

% audio data loading
num_piece = 2;
temp =  audioread("input/9w3r7y_1.wav");
ss{1} = temp(:,1);
[temp, fs] = audioread("input/9w3r7y_2.wav");
ss{2} = temp(:,1);
s = [ss{1}; ss{2}];
num_mic = 1;

% ir data loading
[ir, ~] = audioread("BFTLSPACELIBRARY_ChamberAL.wav");
%[ir, fs] = audioread("logic03sWoodenBooth8000.wav");
%[ir{2}, fs(2)] = audioread("../data/BF TL SPACE LIBRARY/Drumbrella/Drumbrella 5'.R.wav");

%% slra options
useSLRA = 0;
useAmplitude = 1;
avoidNearZeroPolys = 1;
opt.MaxIterations = 100;
opt.epsgrad = 1e-12;
opt.solver = "c"

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
    ir = ir(1:fs);
    L(i) = size(F(ir), 2);

    for j = 1:num_piece
        obs{i, j} = conv(ss{j}, ir);
        X{i, j} = F(obs{i, j});
        N(i, j) = size(X{i, j}, 2);
    end
end
obsCat = [obs{1}; obs{2}];  obsCat = obsCat / max(abs(obsCat));
d = L -1;

%% deconvolution in TF domain
% malloc
Fq = size(X{1}, 1);
for i = 1:2, l(i) = N(1, i) - L(1) + 1; end
RET = zeros(Fq, sum(l));
RET_slra = zeros(Fq, sum(l));
% minimum to 2nd minimum singular value ratio
singValRatio = zeros(Fq, 1);

% frequency-wise processing
parfor f = 1:Fq
    Xf = cell(num_piece,1);
    temp = X{1, 1};
    X1f = temp(f, :)';
    a1 = norm(X1f);
    temp = X{1, 2};
    X2f = temp(f, :)';
    a2 = norm(X2f);

    % (not) structured low-rank approximation
    Syl = [convmtx(X2f, l(1)) convmtx(-X1f, l(2))];
    [~, Sigma, V] = svd(Syl, "econ");
    Sigma = diag(Sigma); singValRatio(f) = Sigma(end) / Sigma(end-1);
    Vtemp = V(:,end);
    v1 = Vtemp(1:l(1)); v1 = a1 * v1 / norm(v1); % scale recovery
    v2 = Vtemp(l(1)+1:end); v2 = a2 * v2 / norm(v2); % scale recovery
    A = [convmtx(v1, L(1)); convmtx(v2, L(1))];     b = [X1f; X2f];
    Hhat = A \ b; % calc. IR (TF domain)
    HhatPhase = Hhat / abs(Hhat);
    [~, idx_HhatPeak] = max(abs(Hhat));
    v1 = v1 / Hhat(idx_HhatPeak); v2 = v2 / Hhat(idx_HhatPeak); % phase normalization
    RET(f,:) = [v1; v2]';

    % agcd by image representation
    averagePower = (a1^2 + a2^2) / sum(N);
    avoidNearZeroPolys = (averagePower > 1);

    if useSLRA % && avoidNearZeroPolys
        if useAmplitude
            A1f = abs(X1f); A2f = abs(X2f);
            [ph, info] = gcd_nls({A1f, A2f}, [], d, opt);
            Hhat = info.Rh(:);
        else
            % どうやら，カーネルが二次元なのでかなり時間がかかる
            [ph, info] = gcd_nls_complex({X1f, X2f}, [], d, opt);
            Hhat = info.hh(:);
        end
        Syl = fullSyl(ph,d);
        [U, Sigma, V] = svd(Syl, "econ");
        Vtemp = V(:,end);
        v1 = Vtemp(1:l(1)); v1 = a1 * v1 / norm(v1); % scale recovery
        v2 = Vtemp(l(1)+1:end); v2 = a2 * v2 / norm(v2); % scale recovery

        [HhatPeak, idx_HhatPeak] = max(abs(Hhat));
        v1 = v1 / Hhat(idx_HhatPeak); v2 = v2 / Hhat(idx_HhatPeak); % phase normalization
        RET_slra(f,:) = [v1; v2]';
        
        if useAmplitude
            % 単純なLRAで計算した際の位相を貼り付ける
            RET_f = RET(f,:);
            RET_f_phase = RET_f / abs(RET_f);
            RET_slra(f,:) = RET_slra(f,:) .* RET_f_phase;
        end
    end
end
%% 
ret = F.pinv(RET);
ret_slra = F.pinv(RET_slra);

F.plotReassign(s);title("dry")
F.plotReassign(obsCat); title("wet")
F.plotReassign(ret); title("lra");
if useSLRA, F.plotReassign(ret_slra); title("slra");end

audiowrite("output/wet.wav", obsCat/max(abs(obsCat)), fs)
audiowrite("output/lra.wav", ret/max(abs(ret)), fs)
audiowrite("output/slra.wav", ret_slra/max(abs(ret_slra)), fs)
% player = audioplayer(obsCat, 8000, 16, 4); play(player);