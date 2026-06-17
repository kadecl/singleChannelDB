F = DGTtool("windowLength", 512, "windowShift", 128);
% speech_data = "/data/p271_006.wav";
path = "input/4077_13754_2sec_11025";
fpaths = dir(path);
speech_data = (path + "/" + fpaths(3).name);

ir_data = "BFTLSPACELIBRARY_ChamberAL_11025.wav";

speech_restore = dereverb(speech_data, ir_data);
speech_restore = F.pinv(speech_restore);

F.plot(speech_restore)
%saveas(gcf, fullfile("/results", "speech_restore.png"));