function rt = reverb_time(x, fs, draw_curve, T)
%REVERB_TIME   Compute Reverberation Time (RT60)
%   rt = reverb_time(x, fs);
%   Compute RT60 from T30 (least square fit on reverberation envelope
%   at -5dB and -35dB).
%
%   rt = reverb_time(x, fs, true); draws figure of reverb time
%
%   rt = reverb_time(x, fs, true, [t1 t2 t3]); draws figure of reverb
%   time with least square fit on reverb envelope from t1 dB to t2 dB
%   and extrapolate the line to t3 dB (by default, t1=-5, t2=-35, t3=-60).
%
%   2008-11-17 MARUI Atsushi

if nargin>3
  T1 = T(1);
  T2 = T(2);
  T3 = T(3);
else
  T1 = -5;
  T2 = -35;
  T3 = -60;
end


re = reverb_envelope(x);
mm = max(10*log10(abs(re)));
t = ([0:length(re)-1]/fs)';
z = 10*log10(abs(re))-mm;
ind = [sum(z >= T1)  sum(z >= T2)];

% linear fit (-5dB ~ -35dB) and find the time of being -60dB
p = polyfit(z(ind(1):ind(2)), t(ind(1):ind(2)), 1);
rt = polyval(p, T3);

if nargin>2 && draw_curve==true
  p = polyfit(t(ind(1):ind(2)), z(ind(1):ind(2)), 1);
  hold on;
  h = plot(t, polyval(p, t), 'k--');
  set(h, 'Color', [.5 .5 .5]);
  title(sprintf('Reverberation Envelope (RT60=%.3fsec)', rt));
  hold off;
end
end

function y = reverb_envelope(x)
% REVERB_ENVELOPE   Calculates reverberation envelope
%    y=reverb_envelope(x) calculates reverberation envelope y from
%    an impulse response x.  This function can be used to see how
%    the reverberation decays in the corresponding room.
%
%    Reference: Ooga, Yamazaki, and Kaneda "Acoustic Systems and
%               Digital Processing for Them" (p.163)
%
% 2007-05-28 MARUI Atsushi

if ~ismatrix(x)
  return;
end

if size(x,1) < size(x,2)
  x = x';
end

y = flipud(cumsum(flipud(x .^ 2)));
end