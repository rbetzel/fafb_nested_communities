%% read in data
% load connectome in matlab
load ../data/connectome_no_threshold.mat

% read numpy community assignments
path2npy = '../data/final_partition_no_threshold.npy';
a = readNPY(path2npy);
for i = 1:size(a,2)
    [~,~,a(:,i)] = unique(a(:,i));
    a(:,i) = fcn_sort_communities(a(:,i));
end
ci = a;

%% calculate coarse-grained community features
filename = '../outputs/community_info.mat';
check = dir(filename);
if isempty(check)

    % generate coarse-grained connectivity
    sc = W.TOT;
    c = dummyvar(ci(:,1));
    d = sum(c)'*sum(c);
    synapse_count = c'*(sc*c);
    synapse_count_norm = synapse_count./d;

    % get downsampled community matrix and optimal ordering
    cidown = zeros(max(ci(:,1)),size(ci,2));
    for i = 1:max(ci(:,1))
        idx = find(ci(:,1)==i,1,'first');
        cidown(i,:) = ci(idx,:);
    end
    dfull = agreement(cidown(:,max(cidown) > 1))/sum(max(cidown) > 1);
    keep = max(cidown) > 3;
    d = agreement(cidown(:,keep))/sum(keep);
    e = 1 - d;
    pd = e(tril(ones(length(e)),-1) > 0);
    ll = linkage(pd','average');
    orderagreement = optimalleaforder(ll,pd);

    order = [];
    tgts = cidown(orderagreement,2);
    while ~isempty(tgts)
        idx = find(cidown(:,2) == tgts(1));
        if length(idx) > 2
            mat = synapse_count(idx,idx);
            pd = pdist(mat + mat','correlation');
            ll = linkage(pd,'average');
            ordercommunity = idx(optimalleaforder(ll,pd));
        end
        order = [order; idx];
        tgts(tgts == tgts(1)) = [];
    end

    % assign colors to communities at each hierarchical level
    [c,cols,ordercols] = fcn_hierarchical_cmap(cidown(:,keep),order);

    % generate a dendrogram
    cikeep = cidown(:,keep);
    cikeepordered = cikeep(order,:);
    treecoor = [];
    for j = 1:size(cikeep,2) - 1
        unq = unique(cikeep(:,j));
        for k = 1:length(unq)
            kdx = find(cikeepordered(:,j) == unq(k));
            par = cikeepordered(kdx,j + 1);
            kdx_plus = find(cikeepordered(:,j + 1) == par(1));
            treecoor = [treecoor; j - 0.5, mean(kdx), j; j + 0.5, mean(kdx_plus), j; nan,nan,j];
        end
        drawnow;
    end
    %%
    save(filename,...
        'synapse_count','synapse_count_norm','ordercommunity',...
        'cidown','keep','dfull','d','e','orderagreement','cikeep',...
        'c','cols','ordercols','order',...
        'treecoor');
    %%
else
    load(filename);
end

%% Figure 2a.

% generate colormap
cmap = [[46 49 146]/255; [249,237,50]/255; 1 0 0];
cmap = interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,256));

% reorder the downsampled community labels
cc = cidown(order,1);
orderneuron = [];
lblsneuron = [];
for i = 1:length(cc)
    idx = find(ci(:,1) == cc(i));
    orderneuron = [orderneuron; idx];
    lblsneuron = [lblsneuron; ones(length(idx),1)*i];
end

% find all edges and their synapse counts
[u,v,wght] = find(W.TOT(orderneuron,orderneuron));
[wght,idxsortwght] = sort(wght,'ascend');
u = u(idxsortwght); v = v(idxsortwght);

% retain all synapses (can vary minsynapse to make a more ``exclusive''
% plot)
minsynapse = 1;
keep = wght >= minsynapse;

% draw the adjacency matrix with synapse counts on log scale
n = length(orderneuron);
f = fcn_rickplot([2,2,4,4]);
scatter(u(keep),n - v(keep) + 1,fcn_sz(log10(wght(keep)),[1.5,10]),log10(wght(keep)),'filled')
axis image;
colormap(cmap);

%% Figure 2d.

% draw the coassignment matrix
f = fcn_rickplot([2,2,4,4]);
imagesc(d(order,order) + eye(length(d)));
colormap(cmap)
colorbar;
axis image off;

%% Figure 2e.

% calculate matching index -- a measure of input/output similarity --
% based on a binary, fine-scale community matrix of synapse counts
% normalized by community size

[mi,mo,mboth] = matching_ind(synapse_count_norm > 0);
mboth = mboth + mboth';

f = fcn_rickplot([2,2,2,2]);
imagesc(mboth(order,order));
colormap(cmap);
axis image;

%% Figure 2f.

% draw euclidean distance between community centroids
load ../data/coordinates.mat
cent = fcn_dummyvar(ci(:,1),[],1)'*coor;
euc = squareform(pdist(cent));

f = fcn_rickplot([2,2,2,2]);
imagesc(euc(order,order));
colormap(cmap);
axis image;

%% Figure 2g.

% draw scatterplot of matching index versus euclidean distance
nbins = 10;
bins = linspace(min(euc(:)),max(euc(:)),nbins + 1);
bins(end) = bins(end) + 1;

mask = zeros(size(euc));
for i = 1:nbins
    mask(euc >= bins(i) & euc < bins(i + 1)) = i;
end
[f,ph] = fcn_boxpts(mboth(~eye(length(euc))),mask(~eye(length(euc))),ones(nbins,3)*0.5,500);
for i = 1:10
    delete(ph.pointshandle{i})
end
mm = rand(size(euc)) < 0.01 & ~eye(length(euc));
sh = scatter(nbins*euc(mm)/max(euc(:)),mboth(mm),'.');
uistack(sh, 'bottom');