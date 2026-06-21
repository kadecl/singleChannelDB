function speech_restore = dereverb(speech_data, ir_data)
F = DGTtool("windowLength", 512, "windowShift", 128);
[speech, Fs] = audioread(speech_data);
[b,a] = cheby2(8, 60, pi * 50 / Fs, 'high');
speech = filter(b, a, speech);

ir = audioread(ir_data);
% ir = ir(1 : Fs*0.37);

Fq = size(F(speech), 1);
N = size(F(speech), 2);
reverb_speech = F(cconv(speech, ir, length(speech)));

IR = F(ir);
M = size(IR, 2);
% M = 50;

rho = 400;
speech_restore = zeros(Fq, N);

parfor f = 1 : Fq
    reverb_speech_bin = reverb_speech(f, :).';
    regul = norm(reverb_speech_bin);
    reverb_speech_bin = reverb_speech_bin./regul;

    speech_restore_bin = admm_algo(reverb_speech_bin, rho, M, N);

    speech_restore_bin = speech_restore_bin .* regul;
    speech_restore(f, :) = speech_restore_bin.';
end
end
