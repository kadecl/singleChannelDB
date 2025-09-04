%%% 20250228 demo.m に，SylvesterSLRAを使う.
%%% ためしに，ちょっとだけ低ランク近似させる

clear all
%% settings
rtdir = "./input2/";
folders = dir(rtdir);
fId = fopen(rtdir + "main20250902.txt", 'w');

% DGT tool
F = DGTtool("windowLength", 512, "windowShift", 128);
fsRe = 11020; fprintf(fId, "fsResample = %d\n", fsRe);
% audio data loading
num_piece = 2;


% フォルダのみを抽出し、'.' や '..' を除外
folderNames = {folders([folders.isdir] & ~ismember({folders.name}, {'.', '..'})).name};

num_mic = 1;

% ir data loading
filenameIR = "BFTLSPACELIBRARY_ChamberAL_11025.wav";
[ir, ~] = audioread(filenameIR);
%[ir, fs] = audioread("logic03sWoodenBooth8000.wav");
%[ir{2}, fs(2)] = audioread("../data/BF TL SPACE LIBRARY/Drumbrella/Drumbrella 5'.R.wav");
fprintf(fId, "IR::  " + filenameIR + "\n\n");

for idx = 1:numel(folderNames)
    try
        % 対象フォルダのパスを構築
        targetFolder = fullfile(rtdir, folderNames{idx});
        fileList = dir(targetFolder);

        % 音声ファイルの格納用セル配列
        ss = {};

        % ファイルごとの処理
        for i = 1:numel(fileList)
            fileName = fileList(i).name;

            % 無視するファイルやフォルダをスキップ
            if fileList(i).isdir || strcmp(fileName, '.DS_Store'), continue; end

            % 音声ファイルの読み込み
            filePath = fullfile(targetFolder, fileName);
            [audioData, fs] = audioread(filePath);

            % ステレオ → モノラル変換（チャンネルを合成）
            if size(audioData, 2) > 1, audioData = sum(audioData, 2); end

            % リサンプリングして保存
            ss{end+1} = resample(audioData, fsRe, fs);
        end

        % 音声ファイル数のチェック
        if numel(ss) ~= 2
            disp("Invalid directory structure: 音声ファイルが2つではありません");
            break;
        end

        s = [ss{1}; ss{2}];

        %% slra options
        useSLRA = 0;
        useAmplitude = 1;
        avoidNearZeroPolys = 1;

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
            ir = ir(1:floor(fs*rt)); ir = ir / norm(ir);
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
        for f = 1:Fq
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

        end
        %%
        ret = cell(1);
        ret{1} = F.pinv(RET1(:,:));  ret{2} = F.pinv(RET2(:,:));
        ret_slra{1} = F.pinv(RET_slra1(:,:));  ret_slra{2} = F.pinv(RET_slra2(:,:));

        % ret_slra = F.pinv(RET_slra);
        len_ss = length(ss{i});

        [SDRpre,SIR,SAR,perm] = bss_eval_sources([obs{1}(1:len_ss); obs{2}(1:len_ss)]', s');
        fprintf("    obs: SDR %2.3f, SIR %2.3f, SAR %2.3f\n", SDRpre, SIR, SAR)

        [SDRpost,SIR,SAR,perm] = bss_eval_sources([ret{1}(1:len_ss); ret{2}(1:len_ss)]', s');
        fprintf("    ret: SDR %2.3f, SIR %2.3f, SAR %2.3f\n", SDRpost, SIR, SAR)

        fprintf(fId, "%s, SDR improvement = %2.3f\n", folderNames{idx}, SDRpost - SDRpre);

        %audiowrite("output/" + sprintf(audiofilename +"_wet_ORACLE.wav", i), obs{i}, fs)
        %audiowrite("output/" + sprintf(audiofilename +"_lra_ORACLE.wav", i), ret{i}, fs)
    catch ME
        disp(['エラー発生: idx = ' num2str(idx) ', メッセージ: ' ME.message]);
        continue;  % 次のループへ
    end
end
fprintf("\n")
fclose(fId);