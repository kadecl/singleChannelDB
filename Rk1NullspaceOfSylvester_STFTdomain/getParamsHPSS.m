function params = getParamsHPSS(F, fs)

%fs
params.fs = fs;

%winlen
params.wlen = F.FFTnum;

%winfunction
params.win = F.win;

%overlaplength
params.overlaplen = F.shift;

%fftlength
params.fftlen = 2^nextpow2(numel(params.win) + 1);

%timefilterlength in samples
timeFilterLength = 0.2;
params.medtfillen = timeFilterLength/((numel(params.win) - params.overlaplen)/params.fs);

%frequencyfilterlength in samples
frequencyFilterLength = 400;
params.medffillen = frequencyFilterLength/(params.fs/params.fftlen);

end