F = DGTtool("windowLength", 512, "windowShift", 128);
speech_data = "/data/p271_006.wav";
ir_data = "/data/IR_sweep_15s_45Hzto22kHz_FS16kHz.v370.wav";

speech_restore = dereverb(speech_data, ir_data);
speech_restore = F.pinv(speech_restore);

F.plot(speech_restore)
saveas(gcf, fullfile("/results", "speech_restore.png"));