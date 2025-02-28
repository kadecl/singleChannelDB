function ph = orthoProjOntoSylvester(S, n, d)
%%% n: degree of polynoials
%%% d: GCD degree

l = n - d + 1; % length of the coefficient vector of quotients
subS{1} = S(:,1:l(1));
subS{2} = S(:,l(1)+1:end);
ph = cell(2,1);
for k = 1:2
    ll = l(3-k);  % n1 : l2 // n2 : l1
    ph_ = zeros(n(k)+1,1);
    for i = 1:ll
        temp = subS{k}(i:i+n(k), i);
        ph_ = ph_ + temp;
    end
    ph{k} = ph_ / ll;
end
