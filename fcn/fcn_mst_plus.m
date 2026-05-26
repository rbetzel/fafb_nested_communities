function [c,amst] = fcn_mst_plus(a,dens)

% amst = graphminspantree(sparse(1 - a),'method','kruskal');
g = graph(double(a + a'));
m = minspantree(g,'method','dense');
amst = m.adjacency;
amst = amst + amst';

n = length(a);
m = n*(n - 1)/2;
tgt = round(m*dens);
mmst = nnz(triu(amst));
add = tgt - mmst;
[u,v,w] = find(a.*triu(a & ~amst,1));
[~,idxsort] = sort(w,'descend');
idx = (v - 1)*n + u;
b = amst;
b(idx(idxsort(1:add))) = 1;
c = double(b | b');