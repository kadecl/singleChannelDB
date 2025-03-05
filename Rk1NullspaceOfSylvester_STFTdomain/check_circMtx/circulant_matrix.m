function C = circulant_matrix(x, m)
    l = length(x);
    if nargin == 1
    m = l;
    end
    C = zeros(l, m);
    for i = 0:m-1
        C(:, i+1) = circshift(x, i);
    end
end
