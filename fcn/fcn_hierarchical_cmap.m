function [c,cols,ordercols] = fcn_hierarchical_cmap(ci,order)
% given a hierarchy of communities, assign distinguishable colors to the
% finest scale and then coloar coarser communities as the mean of color of
% their children

[nr,nc] = size(ci);
cols = distinguishable_colors(nr);
%%
pd = pdist(cols);
ll = linkage(pd,'average');
ordercols = optimalleaforder(ll,pd);
%%
colsordered = cols(ordercols,:);
%%
c = ones([size(ci),3])*0.75;
%%
for i = 1:nc
    unq = unique(nonzeros(ci(:,i)));
    ncomms = length(unq);
    for icomm = 1:ncomms
        idx = ci(order,i) == unq(icomm);
        jdx = ci(:,i) == unq(icomm);
        d = nanmean(colsordered(idx,:),1);
        c(jdx,i,:) = repmat(d,[sum(jdx),1,1]);
    end
end