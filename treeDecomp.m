function [T_parent, T_row, T_col] = treeDecomp(M)
%TREEDECOMP  Returns a tree decomposition of the bipartite adjacency graph
%
%   [Tparent, Trow, Tcol] = treeDecomp(M)
%   Obtain a tree decomposition T = (Tparent, Trow, Tcol) of the bipartite
%   adjacency graph of M, using Approximate Minimum Degree
%   This decomposition is NOT optimal!
%
%   Input:
%       M --> square matrix
%
%   Output:
%       Tparent --> parent pointers of the tree
%       T_row --> cell of row bags
%       T_col --> cell of column bags
%
%   Author:
%       Diego Cifuentes    (diegcif at mit dot edu)

[n,n2] = size(M);
if n~=n2; error('Matrix has to be square');end
nT = 2*n;
T_parent = zeros(1,nT);
T_col = cell(1,nT);
T_row = cell(1,nT);
I = (M~=0);
I = [false(n,n) I; I' false(n,n)];
P = amd(I);
[~,Pinv] = sort(P);
for i=1:nT
    pi = P(i);
    Ni = reshape(find(I(:,pi)),1,[]);
    Ni = Ni(Pinv(Ni)>i);
    if ~isempty(Ni)
        T_parent(i) = min(Pinv(Ni));
    end
    I(Ni,Ni) = true(length(Ni));
    I(pi,:) = 0;
    Ni = union(Ni,pi);
    T_row{i} = Ni(Ni<=n);
    T_col{i} = Ni(Ni>n) - n;
end