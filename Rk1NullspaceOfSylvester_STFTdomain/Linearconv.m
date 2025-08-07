function D =  Linearconv(X,H)
% strict linear convolution in TF domain
% X : 2D
% H : 3D

M=size(X,1);
% N=size(X,2);
D_half=zeros(M/2+1 , size(X,2)+size(H,2)-1);
for m =1:M/2+1
    for q = 1:M
        D_half(m,:) = D_half(m,:) + conv(X(q,:),squeeze(H(m,:,q)));
    end
end
D = [D_half; conj(D_half(end-1:-1:2,:))];
end