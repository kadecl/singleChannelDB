l = [5, 5];
s1 = [1 -2 0 0 4]';
s2 = [1 2 1 3 1]';
h = [1  0 1]; h = h(:);
lh = length(h);

x1 = cconv(s1, h);
x2 = cconv(s2, h);

% check commutativity
norm( cconv(x1, s2) - cconv(x2, s1) )

cx1 = x1;
rx1 = [x1(1) flipud(x1(2:l(2)))'];
%Cx1 = toeplitz(cx1, rx1);
Cx1 = circulant_matrix(x1, length(x2)-1);

e1 = norm(Cx1*s2 - cconv(x1, s2, length(x1)))

cx2 = x2;
rx2 = [x2(1) flipud(x2(2:l(1)))'];
%Cx2 = toeplitz(cx2, rx2);
Cx2 = circulant_matrix(x2, length(x2)-1);
s1ext = buffer(s1, length(x2-1));
e2 = norm(Cx2*s1ext - cconv(x2, s1, length(x2)) )
s = [s1ext; s2ext]; s = s/ norm(s);

commuteness = norm(Cx2*s1 - Cx1*s2)
C = [Cx2 , -Cx1];
C * s
[~,sv,V] = svd(C);
sv = diag(sv); minSingVal =sv(end)
V' * s %%% 零空間が大きいので，その中に埋もれてしまっている