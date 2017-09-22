%   Permanent of a band matrix
%
%   Computes the permanent of a band matrix, using the decomposition 
%   algorithm from
%       "An efficient tree decomposition method for permanents and mixed discriminants"
%       D Cifuentes, P A Parrilo, arxiv.org/abs/1507.03046
%   The complexity is O(n 2^{w1+w2}), where w1, w2 are the lower and upper
%   bandwidths of the matrix.
%
%   Author:
%       Diego Cifuentes    (diegcif at mit dot edu)

% Parameters: lenght n and bandwidths w1, w2
n = 20;
w1 = 5;
w2 = 3;

% Generate random band matrix 
M = randn(n);
for i=1:n
    M(i+w1+1:end,i) = 0;
    M(i,i+w2+1:end) = 0;
end
subplot(1,2,1)
spy(M); xlabel('Sparsity structure')

% Obtain tree decomposition
% Define the tree
nT = 2*n - w1 - w2;
T_parent = zeros(1,nT);
root = ceil(nT/2);
for i=1:(root-1)
    T_parent(i) = i+1;
end
for i=(root+1):nT
    T_parent(i) = i-1;
end

% Define row and column bags
T_col = cell(1,nT);
T_row = cell(1,nT);
for i=1:n-w1-w2
    T_col{w1+2*i-1} = i:i+w1+w2-1;
    T_col{w1+2*i} = i+1:i+w1+w2;
    T_row{w1+2*i-1} = i+w1;
    T_row{w1+2*i} = i+w1;
end
for i=1:w1
    T_col{i} = 1:i+w2;
    T_row{i} = i;
end
for i=1:w2
    T_col{nT-i+1} = n-i-w1+1:n;
    T_row{nT-i+1} = n-i+1;
end
subplot(1,2,2)
treeplot(T_parent); xlabel('Tree decomposition')
pause(.2)

% Compute the permanent using tree decomposition
fprintf('Permanent using tree decomposition\n')
perm = sparsePerm(M,T_parent,T_row,T_col);

% Compute the permanent using Ryser formula
if n <= 20
    fprintf('------------------------------\n')
    fprintf('\nPermanent using Ryser formula\n')
    perm2 = densePerm(M);
    fprintf('Permanent: %0.8g \n', perm2)
end