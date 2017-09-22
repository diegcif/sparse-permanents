function perm = densePerm(M,n)
%DENSEPERM(M) Computes the permanent of a small dense matrix
%
%   perm = densePerm(M) 
%   Computes the permanent of a matrix M using Ryser formula
%   The complexity is O(n 2^n)
%
%   Input:
%       M --> square matrix
%
%   OUTPUT:
%       perm --> permanent of M
%
%   Author:
%       Diego Cifuentes    (diegcif at mit dot edu)

if nargin < 2
    [n,n2] = size(M);
    if n~=n2; error('Matrix has to be square');end
    if ~isfloat(M); error('Input has to be floating point'); end
end

% Use expansion por small values
if n<=4
    if n==0; perm = 1;
    elseif n==1; perm = M;
    elseif n==2; perm = M(1)*M(4)+M(2)*M(3);
    elseif n==3; 
        perm = M(3)*M(5)*M(7)+M(2)*M(6)*M(7)+M(3)*M(4)*M(8)+M(1)*M(6)*M(8)+M(2)*M(4)*M(9)+M(1)*M(5)*M(9);
    else
        perm = M(4)*M(7)*M(10)*M(13)+M(3)*M(8)*M(10)*M(13)+M(4)*M(6)*M(11)*M(13)+M(2)*M(8)*M(11)*M(13)+...
            M(3)*M(6)*M(12)*M(13)+M(2)*M(7)*M(12)*M(13)+M(4)*M(7)*M(9)*M(14)+M(3)*M(8)*M(9)*M(14)+...
            M(4)*M(5)*M(11)*M(14)+M(1)*M(8)*M(11)*M(14)+M(3)*M(5)*M(12)*M(14)+M(1)*M(7)*M(12)*M(14)+...
            M(4)*M(6)*M(9)*M(15)+M(2)*M(8)*M(9)*M(15)+M(4)*M(5)*M(10)*M(15)+M(1)*M(8)*M(10)*M(15)+...
            M(2)*M(5)*M(12)*M(15)+M(1)*M(6)*M(12)*M(15)+M(3)*M(6)*M(9)*M(16)+M(2)*M(7)*M(9)*M(16)+...
            M(3)*M(5)*M(10)*M(16)+M(1)*M(7)*M(10)*M(16)+M(2)*M(5)*M(11)*M(16)+M(1)*M(6)*M(11)*M(16);
    end
    return
end

% Ryser algorithm
perm = 0;
N = 0:pow2(n)-1;
graytable = bitxor(N,bitshift(N,-1));
xold = 0;
sgn = -1;
srow = zeros(n,1);
for x = graytable(2:end) 
    diff = x - xold;
    if diff > 0; 
        [~,j] = log2(diff);
        srow = srow + M(:,j);
    else
        [~,j] = log2(-diff);
        srow = srow - M(:,j);
    end
    perm = perm + sgn * prod(srow,1);
    xold = x;
    sgn = -sgn;
end
perm = perm * (-1)^n;
