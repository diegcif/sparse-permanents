function c = subsetConv(f,g,n)
%SUBSETCONV  Subset convolution of two row vectors
%
%   c = subsetConv(f,g)
%   Performs the subset convolution of two row vectors
%   using the algorithm in
%       Bjorklund, et al. arXiv:cs/0611101
%
%   Input:
%       f, g --> row vectors
%
%   Output:
%       c --> subset convolution of f, g
%
%   Author:
%       Diego Cifuentes    (diegcif at mit dot edu)

if nargin < 3
    [r,N] = size(f);
    [r2,N2] = size(g);
    if N~=N2||r~=1||r2~=1; error('Input must be row vectors');end
    n = log2(N);
    if floor(n)<n; error('Lenght must be a power of 2');end
end

if n > 32; error('Input is too large'); end
N = pow2(n);
wHam = uint32Hamming(0:N-1);
F = rankmobius(f,wHam,n,N);
G = rankmobius(g,wHam,n,N);
C = convCols(F,G,n,N);
c = invrankmobius(C,wHam,n,N);

function w = uint32Hamming( x )
% Hamming weight

x = uint32(x);
x = bitand(bitshift(x, -1), uint32(1431655765)) + bitand(x, uint32(1431655765));
x = bitand(bitshift(x, -2), uint32(858993459)) + bitand(x, uint32(858993459));
x = bitand(bitshift(x, -4), uint32(252645135)) + bitand(x, uint32(252645135));
x = bitand(bitshift(x, -8), uint32(16711935)) + bitand(x, uint32(16711935));
x = bitand(bitshift(x, -16), uint32(65535)) + bitand(x, uint32(65535));
w = double(x);

function F = rankmobius(f,wHam,n,N)
% Ranked Mobius transform

F = full(sparse(wHam(1:N)+1,1:N,f));
for i=2:n+1
    Fold = F;
    j = i-1; 
    J0 = repmat(0:pow2(j-1)-1,[1,pow2(n-j)]);
    J1 = reshape(repmat(pow2(j)*(0:pow2(n-j)-1),[pow2(j-1),1]),1,[]);
    J = 1 + pow2(j-1) + J0 + J1; %indices have j
    Jc = 1 + J0 + J1; %indices dont have j
    F(:,Jc) = Fold(:,Jc);
    F(:,J) = Fold(:,J) + Fold(:,Jc);
end

function f = invrankmobius(F,wHam,n,N)
% Inverse ranked Mobius transform

for i=2:n+1
    Fold = F;
    j = i-1; 
    J0 = repmat(0:pow2(j-1)-1,[1,pow2(n-j)]);
    J1 = reshape(repmat(pow2(j)*(0:pow2(n-j)-1),[pow2(j-1),1]),1,[]);
    J = 1 + pow2(j-1) + J0 + J1; %indices have j
    Jc = 1 + J0 + J1; %indices dont have j
    F(:,Jc) = Fold(:,Jc);
    F(:,J) = Fold(:,J) -Fold(:,Jc);
end
I = sub2ind(size(F), wHam(1:N)+1, 1:N);
f = F(I);

function C = convCols(A,B,n,N)
% Convolution of the columns of A, B

if n < 15
    C = zeros(n+1,N);
    for k=0:n
        for j=0:k
            C(k+1,:) = C(k+1,:) + A(j+1,:).*B(k-j+1,:);
        end
    end
else
    Z = zeros(n+1,N);
    C = ifft(fft([A;Z]).*fft([B;Z]));
    C = C(1:n+1,:);
end