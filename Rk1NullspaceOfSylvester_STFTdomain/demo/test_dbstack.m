A = randn(4,4);

st = dbstack('-completenames');
[folderPath, fileName, ext] = fileparts(st(1).file);
resultFileName = strcat(fileName,  '_result')
saveFilePath = fullfile(folderPath, [resultFileName, '.mat']);
data = rand(100, 1);
save(saveFilePath, 'data');
fprintf('結果を %s に保存しました。\n', saveFilePath);