% Configuration
L = 2^10;
x = randn(L,1); %signal
M = 128; %window length
a = M/4; %window shift
c = lcm(a,M);
L = ceil(length(x)/c)*c;

b = L/M;
N = L/a;

fs=M;
r = -1 + 2*rand(fs,1) ;
r(1) = 1;
t=(0:1/fs:1-1/fs)';
h=r.*exp(-10*t); % FIR
h=padarray(h,M,0,"both"); % zero padding
h=buffer(h,ceil(length(h)/c)*c);
K=length(h); % FIR length

F = DGTtool('windowLength',M,'windowShift',a,'windowName','h');
F.makeWindowTight
X = F(x);
Xfull = specHalf2Full(X);
norm(x-F.pinv(X))/norm(x) %reconstruction check
dorig = cconv(x,h,length(x));  %reference in time domain
Dorig = F(dorig); %reference in T-F domain

% I made 2 patterns of computing H3D
% 1. FFT and Hadamard Product Algorithm (See Algorithm 1)
C=ComputeH_Conv(F, length(h), M, length(h)/M);
H3D=ComputeH(C, h, a, M, length(h)/M);


% 2. Overlap-add or Factorization algorithm (See Chapter 4 & 5)
FforJ = DGTtoolOLA('windowLength',M,'windowShift',1,'windowVector',F.dualWin);
FforH = DGTtoolOLA('windowLength',M,'windowShift',a,'windowVector',F.win);


% compute with the factorization algorithm if executed.
% compute with the Overlap-add algorithm if commented-out.
FforH.win =buffer(FforH.win,K);

J=FforJ(h);   % computing eq(17) in the H3DFactorization.pdf
H=FforH(permute(J,[2 1]));  % computing eq(16) in the H3DFactorization.pdf


norm(H3D-H,"fro")/norm(H3D,"fro") % check the same result


REJ=(permute(FforH.pinv(H),[2 1]));
norm(J-REJ,"fro")/norm(J,"fro") % check the reconstruction of J
reh=FforJ.pinv(J);
norm(reh-h)/norm(h) % check the reconstruction of h


D=Strictconv(Xfull,H3D); % convolution in T-F domain
d=F.pinv(D);
difTFdomain = norm(D-Dorig,"fro")/norm(Dorig,"fro")
difTdomain  = norm(d-dorig)/norm(dorig)
norm_T(M)   = norm(D-Dorig,"fro")/norm(Dorig,"fro");
norm_TF(M)  = norm(d-dorig)/norm(dorig);

%% Inverse 
Jre=permute(FforH.pinv(H),[2 1]);
hre=FforJ.pinv(Jre);

norm(J-Jre,"fro")/norm(J,"fro")
norm(h-hre)/norm(h)


hrere = InvComputeH(C,F,H3D,a,M,length(h)/M);
norm(h-hrere)/norm(h)
figure,plot(h)
figure,plot(hrere)

%% eq(5) test: changing the amount of summation of q
% (the convolution is strict when q=0:M-1)
D=StrictconvApp(Xfull,H3D);
d=F.pinv(D);
norm_TF(1)  = norm(D(:)-Dorig(:))/norm(Dorig(:));
norm_T(1)   = norm(d-dorig)/norm(dorig);
for i =3:2:M
    D=StrictconvApps(Xfull,H3D,i);
    d=F.pinv(D);
    norm_TF(i)  = norm(D(:)-Dorig(:))/norm(Dorig(:));
    norm_T(i)   = norm(d-dorig)/norm(dorig);
end


semilogy([3:2:M M],norm_TF([3:2:M M]))
hold on
semilogy([3:2:M M],norm_T([3:2:M M]))
hold off
%%

function Z = specHalf2Full(X)
Z = [X; conj(flipud(X(2:end-1,:)))];
end

function X = specFull2Half(Z)
X = Z(1:floor(size(Z,1)/2)+1,:);
end



function D=  Strictconv(Xfull,H)
M=size(Xfull,1);
N=size(Xfull,2);
D=zeros(M/2+1,N);
for m=0:M/2
    for q=0:M-1
        D(m+1,:) =  D(m+1,:) + cconv(Xfull(q+1,:),H(m+1,:,q+1),N);
    end
end
end

function D=  StrictconvApp(Xfull,H)
% This code only computes subband convolution in the T-F domain.
M=size(Xfull,1);
N=size(Xfull,2);
D=zeros(M/2+1,N);
for m=0:M/2
    D(m+1,:) =  D(m+1,:) + cconv(Xfull(m+1,:),H(m+1,:,m+1),N);
end
end

function D=  StrictconvApps(Xfull,H,s)
% This code computes crossband convolution in the T-F domain.
% But this is the approximated computation.
M=size(Xfull,1);
N=size(Xfull,2);
D=zeros(M/2+1,N);
for m=0:M/2
    D(m+1,:) =  D(m+1,:) + cconv(Xfull(m+1,:),H(m+1,:,m+1),N);
    for q=1:(s-1)/2
        D(m+1,:) =  D(m+1,:) + cconv(Xfull(mod(m-q,M)+1,:),H(m+1,:,mod(m-q,M)+1),N);
        D(m+1,:) =  D(m+1,:) + cconv(Xfull(mod(m+q,M)+1,:),H(m+1,:,mod(m+q,M)+1),N);
    end
end
end
