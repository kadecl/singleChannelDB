%%% STFT domain では， 周波数領域への変換で波形のcyclic性を仮定しているので，
%%% 少なくともframeの内部での畳み込みは巡回畳み込みでモデル化されるべきである．
%%% 一方，frame間での畳み込みは巡回性を仮定してしまうのはどうなんだろう..？
%%% frame内部での畳み込みから時間周波数領域上での厳密な畳み込みを導出すると，frame間の畳み込みもcyclicになってしまう．
%%% はじっこに影響出るだけなのでしょうがないよね〜〜って感じ?
%%% frame内が線形畳み込みになってると仮定した時のTF domainげんみつ畳み込みはどうなっているだろうか

l = [5, 5];
s1 = [1 -2 0 0 4]';
s2 = [1 2 1 3 1]';
h = [1  0 1]; h = buffer(h(:), length(s1));
lh = length(h);

x1 = cconv(s1, h);
x2 = cconv(s2, h);

% check commutativity
norm( cconv(x1, s2) - cconv(x2, s1) )

s2ext = s2;%buffer(s2, length(x2)-1);
Cx1 = circulant_matrix(x1, length(s2ext));
e1 = norm(Cx1*s2ext - cconv(x1, s2, length(x1)))

Cx2 = circulant_matrix(x2, length(s2ext));
s1ext = s1; %buffer(s1, length(x1)-1);
e2 = norm(Cx2*s1ext - cconv(x2, s1, length(x2)) )

commuteness = norm(Cx2*s1ext - Cx1*s2ext)
C = [Cx2 , -Cx1];
s = [s1ext; s2ext]; s = s/ norm(s);
resid = norm(C * s)

[~,sv,V] = svd(C);
sv = diag(sv)';
sim = s'*V; %%% 零空間が大きいので，その中に埋もれてしまっている
sizeDiff = max(size(C)) - min(size(C)) 
% Cを縦長にするには，観測をゼロうめした上で..でめんどくさい．
% 今回の設定なら，はみ出た分として計算できる
s' * V(:,end)
