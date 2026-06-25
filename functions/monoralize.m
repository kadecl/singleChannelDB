function s = monoralize(s)
[m,n] = size(s);
if n > m
    s = s';
end
s = s(:, 1);
end