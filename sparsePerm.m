function [perm,PTable] = sparsePerm(M,T_parent,T_row,T_col)
%SPARSEPERM Permanent of a sparse matrix
%
%   [perm,PTable] = sparsePerm(M)
%   Computes the permanent of a sparse matrix M using the decomposition 
%   algorithm from
%       "An efficient tree decomposition method for permanents and mixed discriminants"
%       D Cifuentes, P A Parrilo, arxiv.org/abs/1507.03046
%   The complexity is O(n 2^k), where k is the width of the decomposition.
%   The tree decomposition is generated using Approximate Minimum Degree
%
%   [perm,PTable] = sparsePerm(M,Tparent,Trow,Tcol)
%   Computes the permanent using a given tree decomposition 
%   T = (Tparent,Trow,Tcol) of the bipartite adjacency graph of M
%
%   Input:
%       M --> square matrix
%       T_parent --> parent pointers of the tree (OPTIONAL)
%       T_row --> cell of row bags (OPTIONAL)
%       T_col --> cell of column bags (OPTIONAL)
%
%   Output:
%       perm --> permanent of M
%       PTable --> dynamic programming table computed
%
%   Examples:
%       exmp_sparsematrix.m --> permanent of a random sparse matrix
%       exmp_bandmatrix.m --> permanent of a band matrix
%
%   Author:
%       Diego Cifuentes    (diegcif at mit dot edu)

if nargin==1; [T_parent,T_row,T_col] = treeDecomp(M); end
check_input(M,T_parent,T_row,T_col);
[T_parent,T_row,T_col,order] = preprocess_tree(T_parent,T_row,T_col);
K = cellfun(@(c) length(c),T_row) + cellfun(@(c) length(c),T_col);
n = length(M);
nT = length(T_parent);
fprintf('Process tree decomposition\nNumber of nodes: %d\nMaximum width: %d\n',nT,max(K))

PTable = cell(1,nT);
for j=1:nT
    i = order(j);
    fprintf('Node %d/%d, width %d\n',j,nT,K(i))
    Ai = T_row{i};
    Xi = T_col{i};
    child_i = find(T_parent == i);
    [Ii,Ic] = find_indicators(n,T_row,T_col,i,child_i);
    permsXi = subPerms(M,Ai,Xi);
    if isempty(child_i)
        PTable{i} = permsXi;
    else
        Permschild = PTable(child_i);
        PTable{i} = evalRecursion(permsXi,Permschild,Ii,Ic);
    end
end
perm = full(PTable{i}(1));
fprintf('Completed\nPermanent: %0.8g \n', perm)

function check_input(M,T_parent,T_row,T_col)
%   Basic verification of the input

[n1,n2] = size(M);
if n1~=n2; error('Matrix has to be square');end
nT = length(T_parent);
if length(T_row)~=nT || length(T_col)~=nT
    error('Lenghts of T_row, T_col are different');
end
if ~any(T_parent==0) || any(T_parent>nT) || any(T_parent<0)
    error('Bad vector of parent pointers');end
for i=1:nT
    Ai = T_row{i};
    Xi = T_col{i};
    if ~isempty(Ai) && (any(Ai>n1) || any(Ai<1))
        error('Indices of T_row exceed matrix dimensions');end
    if ~isempty(Xi) && (any(Xi>n2) || any(Xi< 1))
        error('Indices of T_col exceed matrix dimensions');end
    M(Ai,Xi) = 0;
end
if nnz(M)>0; error('Indices of T_row, T_col do not cover the whole matrix');end

function [T_parent,T_row,T_col,order] = preprocess_tree(T_parent,T_row,T_col)
%   Makes the tree connected and obtains a topological ordering

nT = length(T_parent);
root = (T_parent==0);
if sum(root) > 1
    nT = nT + 1;
    T_parent(nT) = 0;
    T_parent(root) = nT;
    T_row{nT} = [];
    T_col{nT} = [];
end
order = zeros(1,nT);
remparent = T_parent;
count = 0;
while(count<nT)
    leavs = setdiff((1:nT),remparent);
    L = length(leavs);
    order(count+1:count+L) = leavs;
    count = count + L;
    remparent(leavs) = leavs;
end

function [Ii,Ic] = find_indicators(n,T_row,T_col,i,child_i)
%   Obtains indicator vectors of node i and its children

Ii = set2indic(n,T_row{i},T_col{i});
nc = length(child_i);
Ic = false(nc,2*n);
for j = 1:nc
    c = child_i(j);
    Ac = T_row{c};
    Xc = T_col{c};
    Ic(j,:) = set2indic(n,Ac,Xc);
end

function I_S = set2indic(n,SA,SX)
%   Obtains indicator vector of a pair (SA,SX)

I_A = sparse(ones(1,length(SA)),SA,true,1,n);
I_X = sparse(ones(1,length(SX)),SX,true,1,n);
I_S = [I_A, I_X];

function perms = subPerms(M,Ai,Xi)
%   Computes the permanents of some submatrices of M. The matrices considered 
%   have row set contained in Ai, and column set contained in Xi.

n = length(M);
Isubsets = list_partitionsXiAi(n,Ai,Xi,0);
m = size(Isubsets,1);
perms = [sparse(m,1), Isubsets];
for s=1:m
    k = nnz(Isubsets(s,:))/2;
    subM = M(Isubsets(s,1:n), Isubsets(s,n+1:2*n));
    perms(s,1) = densePerm(subM,k);
end

function Isubsets = list_partitionsXiAi(n,Ai,Xi,offset)
%   Find all pairs (D,Y) such that:
%   (1)  D is subset of Ai
%   (2)  Y is subset of Xi 
%   (3)  |Y|-|D|=offset

ncols = length(Xi);
nrows = length(Ai);
diff_row_col = offset;
indic_v = [];
indic_x = [];
count = 0;
for sc = ncols:-1:0
    sr = sc - diff_row_col;
    if sr > nrows || sr < 0; continue;  end
    colsubsets = nchoosek_subsets(Xi,sc);
    rowsubsets = nchoosek_subsets(Ai,sr);
    nc = size(colsubsets,1);
    nr = size(rowsubsets,1);
    indv = zeros(nc*nr,sr+sc);
    indx = zeros(nc*nr,sr+sc);
    for c=1:nc
        indx(1+(c-1)*nr: c*nr, :) = repmat((c-1)*nr + (1:nr).',[1,sr+sc]);
        indv(1+(c-1)*nr: c*nr, :) = [rowsubsets repmat(colsubsets(c,:)+n,[nr,1])];
    end
    indic_v = [indic_v; indv(:)];
    indic_x = [indic_x; count + indx(:)];
    count = count + nc*nr;
end
Isubsets = sparse(indic_x',indic_v',true,count,2*n);

function subs = nchoosek_subsets(S,k)
%   Lists all subsets of S with k elements

l = length(S);
if (k==0); subs = zeros(1,k);
elseif (l==0); subs = zeros(1,k);
elseif (l==1 && k==1); subs = reshape(S,1,[]);
elseif (l==1); subs = zeros(0,k);
else subs = nchoosek(S,k);
end

function P = evalRecursion(permsXi,Perms,Ii,Ic)
%   Evaluates the recursion formula of the dynamic program.
%   For a pair (D,Y), the recursion is the following:
%   P(D,Y) = \sum_{D,Y} P(D_s,Y_s) 
%                \prod_c (-1)^{|D_cs|} P(D_cs,Y_cs) P(D_cc,Y_cc)

n = length(Ii)/2;
nc = size(Ic,1);
P = permsXi;
Pv = bits2vec(P,Ii,false);
for j = 1:nc
    Icj = Ic(j,:);
    I_cs = Icj & Ii;
    ind_rows = search_pattern(I_cs,permsXi(:,2:end), Ii, true);
    Perms_c = permsXi(ind_rows,:);
    Pc = bits2vec(Perms_c,Ii,true);
    Pv = subsetConv(Pv,Pc,nnz(Ii));
end
for j = 1:nc
    Icj = Ic(j,:);
    ind_rows = search_pattern(true(1,2*n),Perms{j}(:,2:end), Icj&~Ii, false);
    Perms_c = Perms{j}(ind_rows,:);
    Pc = bits2vec(Perms_c,Ii,false);
    Pv = subsetConv(Pv,Pc,nnz(Ii));
end
P = vec2bits(Pv,Ii,n);


function Pv = bits2vec(Pb,Ii,hassign)
%   Converts a list indexed by subsets to a vector

if ~any(Ii); Pv=Pb(1,1); return; end
k = nnz(Ii);
So = Pb(:,2:end);
x = 1 + bin2dec(num2str(So(:,Ii)));
y = Pb(:,1);
if hassign; y = y.*((-1).^(sum(So,2)/2)); end
Pv = sparse(ones(1,size(Pb,1)),x,y,1,2^k);

function Pb = vec2bits(Pv,Ii,n)
%   Converts a vector to a list indexed by subsets

k = nnz(Ii);
Pb = zeros(nnz(Pv),2*n+1);
if ~nnz(Pv); Pb = sparse(1,2*n+1);end
ind = find(Pv);
for i = 1:length(ind)
    ii = length(ind)-i+1;
    Pb(ii,1) = Pv(ind(i));
    Is = zeros(1,2*n);
    Is(Ii) = bitget(ind(i)-1,(k:-1:1));
    Pb(ii,2:end) = Is;
end
Pb = sparse(Pb);

function rows = search_pattern(pattern,table, Icols, subpattern)
%   Searches for a pattern (or subpattern)

if subpattern
    zerocols = ~pattern; 
    Icols = Icols & zerocols;
end
pattern = pattern(:,Icols);
table = table(:,Icols);
rows = all(bsxfun(@eq,pattern,table),2);
