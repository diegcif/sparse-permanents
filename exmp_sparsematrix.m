%   Permanent of a sparse matrix
%
%   Computes the permanent of a random sparse matrix, using the decomposition 
%   algorithm from
%       "An efficient tree decomposition method for permanents and mixed discriminants"
%       D Cifuentes, P A Parrilo arxiv.org/abs/1507.03046
%   The complexity is O(n 2^k),  where k is the width of the decomposition.
%
%   Author:
%       Diego Cifuentes    (diegcif at mit dot edu)

% Parameters: lenght n and sparsity s
n = 20;
s = 3/n;

% Generate random sparse matrix
M = 5*sprand(n,n,s,1/(100*n));
subplot(1,2,1)
spy(M); xlabel('Sparsity structure')

% Obtain a tree decomposition using Approximate Minimum Degree
% This decomposition is NOT optimal!
[T_parent, T_row, T_col] = treeDecomp(M);
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