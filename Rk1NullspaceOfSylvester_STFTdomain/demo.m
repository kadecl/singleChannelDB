% load
F = DGTtool("windowLength", 512, "windowShift", 128);

[s, fs_s] = audioread("../data/brahms_vc_8000Hz.wav");
num_piece = 2;
ss{2} = s(1:round(length(s)/2));
ss{1} = s(round(length(s)/2)+1:end);
num_mic = 1;
[ir{2}, fs(2)] = audioread("../data/BF TL SPACE LIBRARY/Chambers/Chamber A -> B.L.wav");
[ir{1}, fs(1)] = audioread("../data/BF TL SPACE LIBRARY/Drumbrella/Drumbrella 5'.R.wav");

% preprocessing
ir_resmpl = cell(num_mic,1);
rt = zeros(num_mic, 1);
obs = cell(num_mic,num_piece);
X = cell(num_mic, num_piece);
N = zeros(num_mic, num_piece);

% resample and shortening
for i = 1:num_mic
    ir_temp = resample(ir{i}, fs_s, fs(i));

    % calcurate the reververation time RT60
    rt(i) = reverb_time(ir_temp, fs_s);
    ir_resmpl{i} = ir_temp( 1 : ceil(fs_s * rt(i)) );
    L(i) = size(F(ir_resmpl{i}), 2);

    for j = 1:num_piece
        obs{i, j} = conv(ss{j}, ir_resmpl{i});
        X{i, j} = F(obs{i, j});
        N(i, j) = size(X{i, j}, 2);
    end
end

Fq = size(X{1}, 1);
l1 = N(1, 1) - L(1) + 1;
l2 = N(1, 2) - L(1) + 1;
RET = zeros(Fq, l1 + l2);
for f = 1:Fq
    Xf = cell(num_piece,1);
    temp = X{1, 1};
    X1f = temp(f, :)';
    temp = X{1, 2};
    X2f = temp(f, :)';
    Syl = [convmtx(X1f, l2) convmtx(X2f, l1)];
    [U, Sigma, V] = svd(Syl, "econ");
    Vtemp = V(:,end);
    v2 = Vtemp(1:l2);
    v1 = Vtemp(l2+1:end);
    A = [convmtx(v1, L(1)); convmtx(v2, L(1))];
    b = [X1f; X2f];
    Hhat = A \ b;
    PhaseNormalizationCoeff = Hhat(1) / abs(Hhat(1));
    RET(f,:) = Vtemp;% * PhaseNormalizationCoeff; % phase normalization. c.f. yatabe2024icassp
end
%% 

ret = F.pinv(RET);