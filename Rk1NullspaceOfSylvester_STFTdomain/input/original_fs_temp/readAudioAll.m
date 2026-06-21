d = dir('./*.wav');
for i = 1:length(d)
    fname = d(i).name
    [s,fs] = audioread(fname);
    s_re = resample(s, 11025, fs);

    tok = regexp(fname, '^(.*?)(?=_bip)', 'once', 'match', 'ignorecase')
    if isempty(tok)
        continue; % "_bip" が含まれないファイルはスキップ
    end
    dirName = tok            % 取得した文字列
    % 必要なら拡張子や末尾の不要文字を取り除く（例: 空白削除）
    dirName = strtrim(dirName);

    % ディレクトリパス（同じフォルダ内に作る場合）
    targetDir = "../inst/" + dirName;

    % 存在しなければ作成
    if ~exist(targetDir, 'dir')
        mkdir(targetDir);
    end

    audiowrite(targetDir + "/" + fname, s_re, 11025);
end