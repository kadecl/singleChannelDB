function speech_restore = dereverb(X, Lir)
% F = DGTtool("windowLength", 512, "windowShift", 128);
%[speech, Fs] = audioread(speech_data);
%[b,a] = cheby2(8, 60, pi * 50 / Fs, 'high');
%speech = filter(b, a, speech);

% ir = audioread(ir_data);
% ir = ir(1 : Fs*0.37);

[Fq, N] = size(X);
% N = size(X, 2);
% reverb_speech = F(cconv(speech, ir, length(speech)));

M = Lir;
% M = 50;

rho = 400;
speech_restore = zeros(Fq, N);

parfor f = 1 : Fq
    %reverb_speech_bin = reverb_speech(f, :).';
    reverb_speech_bin = X(f, :).';
    regul = norm(reverb_speech_bin);
    reverb_speech_bin = reverb_speech_bin./regul;

    speech_restore_bin = admm_algo(reverb_speech_bin, rho, M, N);

    speech_restore_bin = speech_restore_bin .* regul;
    speech_restore(f, :) = speech_restore_bin.';
end
end
