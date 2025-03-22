function ret = weightedMeanWithZeroPadding(s, x, w)

s = s(:); x = x(:);
ls = length(s); lx = length(x);
if ls > lx
    s = s(1:lx);
else
    x = x(1:ls);
end

if nargin == 3
    w = w(1:min(ls, lx));
    ret = sum( w.* (s-x) / sum(w));
else
    ret = mean(s-x);
end