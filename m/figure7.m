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

% load annotations
load ../data/annotations.mat

%%
load ../data/coordinates.mat
load ../data/annotations.mat
load ../data/visual_neuron_data.mat
load ../data/spatial.mat
%%
for i = 1:size(a,2)
    [~,~,a(:,i)] = unique(a(:,i));
    a(:,i) = fcn_sort_communities(a(:,i));
end
ci = a;
n = size(ci,1);
%% define some useful quantities and load some useful data
load ../outputs/community_info.mat
cmap = [[46 49 146]/255; [249,237,50]/255; 1 0 0];
cmap = interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,256));
tgts = fieldnames(labels);
%%
addpath ../fcn

%%

lvl = 1;
filename = sprintf('../outputs/community_features.lvl=%i.mat',lvl);
load(filename);

ncoarse = max(ci(:,lvl));

features.spatial = rmfield(features.spatial,'soma_centroid');
features.spatial = rmfield(features.spatial,'arbour_centroid');

fn = fieldnames(features);
fn(contains(fn,'neurotransmitter')) = [];

X = [];
N = [];
L = [];
for j = 1:length(fn)
    y = features.(fn{j});
    fnsub = fieldnames(y);
    for k = 1:length(fnsub)
        x = [y.(fnsub{k})];
        if contains(fn{j},'neurotransmitter')
            m = length(x)/length(synapse_count_norm);
            x = reshape(x,[m,ncoarse]);
            % px = bsxfun(@rdivide,x,sum(x,1));
            % entx = -nansum(px.*log2(px),1)/log2(length(px));
            % x = entx;
        end
        if contains(fn{j},'sbm')
            m = length(x)/length(synapse_count_norm);
            x = reshape(x,[m,ncoarse])';
            x = nanmean(x,2);
        end
        [nr,nc] = size(x);
        if nc == length(synapse_count_norm)
            x = x';
        end
        X = [X,x];
        if contains(fn{j},'neurotransmitter')
            N = [N; strcat(ntnames,'_',repmat(fnsub(k),size(x,2),1))];
        else
            N = [N; repmat(fnsub(k),size(x,2),1)];
        end

        L = [L; j*ones(size(x,2),1)];
    end
end
X = [X,assort_data];
L = [L; ones(size(assort_data,2),1)*(j + 1)];

% regs = contains(N,'size');
% Xr = X;
% for i = 1:size(X,2)
%     [~,~,Xr(:,i)] = regress(X(:,i),...
%         [ones(ncoarse,1),...
%         X(:,regs),X(:,regs).^0.5,X(:,regs).^2]);
% end
Xr = X;
z = bsxfun(@rdivide,bsxfun(@minus,Xr,nanmean(Xr)),nanstd(Xr));
N = [N; {'assortativity'}; {'disassortativity'}; {'core'}; {'periphery'}; {'diversity'}];

drop = isnan(sum(full(z),2));

[coeff,score,latent] = pca(full(z));
latent = latent/sum(latent);
clatent = cumsum(latent);
varexp = 0.85;
feat_lodim = score(~drop,1:find(clatent >= varexp,1,'first'));

%%

maxIter = 1000;
kmax = 250;
options = statset('maxiter',10000);
lims = [2,kmax];
niter0 = 25;
inc = 2;
nstages = 4;
niter = 100;
cutoff = 1;
maxinternal = 100000;
stopinternal = round(maxinternal*0.05);
niterconsensus = 100;

ndims = size(feat_lodim,2);

gmmdir = fullfile('..','outputs',sprintf('gmm_varexp=%.4f',varexp));

% UNCOMMENT THIS IF YOU WANT TO FIT GMMs YOURSELF
% % if ~exist(gmmdir,'dir'); mkdir(gmmdir); end
% % %%
% % 
% % aic = nan(niter,kmax);
% % bic = aic;
% % % for k = 2:kmax
% % save(fullfile(gmmdir,'data_no_nts.mat'),'feat_lodim','drop','ndims','maxIter','varexp');
% % for k = 2:kmax
% %     mdl = cell(niter,1);
% %     gmmfilename = fullfile(gmmdir,sprintf('k=%1.4i.data_no_nts.mat',k));
% %     check = dir(gmmfilename);
% %     if isempty(check)
% %         pi_k = zeros(k,niter,'single');
% %         mu_k = zeros(k,ndims,niter,'single');
% %         Sigma_k = mu_k;
% %         logL = cell(niter,1);
% %         clu = zeros(sum(~drop),niter,'single');
% %         resp = zeros(sum(~drop),k,niter,'single');
% %         AIC = zeros(niter,1,'single');
% %         BIC = AIC;
% %         parfor iter = 1:niter
% %             disp([k,iter]);
% %             [~, warnid] = lastwarn;
% %             warning('off', warnid);
% %             [pi_k(:,iter), mu_k(:,:,iter), Sigma_k(:,:,iter), logL{iter}, clu(:,iter), resp(:,:,iter), AIC(iter), BIC(iter)] = fcn_gmm_fit(feat_lodim,k,maxIter);
% %         end
% %         save(gmmfilename,'pi_k','mu_k','Sigma_k','logL','clu','resp','AIC','BIC');
% %     else
% %         load(gmmfilename,'BIC','AIC')
% %     end
% %     aic(:,k) = AIC;
% %     bic(:,k) = BIC;
% % 
% % end
% % 
% % aic_m = nanmean(aic,1);
% % bic_m = nanmean(bic,1);
% % 
% % aic_s = nanstd(aic,[],1);
% % bic_s = nanstd(bic,[],1);
% % 
% % [~,aic_min] = min(aic_m);
% % [~,bic_min] = min(bic_m);
% % 
% % aic_prct = prctile(aic(~isnan(aic)),cutoff);
% % bic_prct = prctile(bic(~isnan(bic)),cutoff);
% % 
% % aic_count = sum(aic <= aic_prct,1);
% % bic_count = sum(bic <= bic_prct,1);
% % 
% % [r,c] = find(bic <= bic_prct);
% % ind = sub2ind(size(bic),r,c);

consensusfilename = fullfile(gmmdir,sprintf('consensus.data_no_nts.bic=%0.4f.mat',cutoff));
check = dir(consensusfilename);
if isempty(check)

    unq = unique(c);
    clu_good = zeros(sum(~drop),length(r));
    count = 0;
    for i = 1:length(unq)
        k = unq(i);
        load(fullfile(gmmdir,sprintf('k=%1.4i.data_no_nts.mat',k)));
        idx = c == unq(i);
        clu_keep = clu(:,r(idx));
        idx = (count + 1):(count + sum(idx));
        count = idx(end);
        clu_good(:,idx) = clu_keep;
    end

    clu_consensus = zeros(size(clu_good,1),niterconsensus);
    fitness_consensus = cell(niterconsensus,1);
    parfor iter = 1:niterconsensus
        disp([v,iter]);
        [clu_consensus(:,iter),~,fitness_consensus{iter}] = fcn_consensus_communities_greedy(clu_good,maxinternal,stopinternal);
    end
    save(consensusfilename,'clu_consensus','fitness_consensus','bic_prct','bic_min','bic_count','r','c','unq','clu_good');
else
    load(consensusfilename)
end

mx = cell2mat(cellfun(@max,fitness_consensus,'uniformoutput',false));
[~,idxbest] = max(mx(:,2));

[~,~,clu_relabel] = unique(clu_consensus(:,idxbest));

clu_ind = dummyvar(clu_good);
d = clu_ind*clu_ind';
d = d/size(clu_good,2);
% d = agreement(clu_good)/size(clu_good,2);
% d = d + eye(length(d));

rho = corr(z(~drop,:)');

cent = z(~drop,:)'*fcn_dummyvar(clu_relabel,[],1);
pd = pdist(cent','correlation');
ll = linkage(pd,'average');
order_cent = optimalleaforder(ll,pd);
cols_cent = distinguishable_colors(length(order_cent));
pd_cols = pdist(cols_cent);
ll_cols = linkage(pd_cols,'average');
order_cols = optimalleaforder(ll_cols,pd_cols);
cols_cent = cols_cent(order_cols,:);


clu_relabel_new = clu_relabel;
for i = 1:length(order_cent)
    clu_relabel_new(clu_relabel == order_cent(i)) = i;
end
clu_relabel = clu_relabel_new;

clubest = zeros(ncoarse,1);
clubest(~drop) = clu_relabel;

% clubest_aggregated(:,v) = clubest;

clu_relabel_manual = clu_relabel;
clu_relabel_manual(clu_relabel == 1) = 17;
clu_relabel_manual(clu_relabel == 51) = 50;

%%

tsne_name = fullfile(gmmdir,'tsne_coordinates.no_nts.mat');
check = dir(tsne_name);
if isempty(check)
    coor_tsne = tsne(feat_lodim,'perplexity',20,'Exaggeration',5);
    save(tsne_name,'coor_tsne')
else
    load(tsne_name)
end
%% Fig 7a

% cols_cent = distinguishable_colors(length(order_cent));

num_nodes = hist(ci(:,1),1:ncoarse);
num_nodes = num_nodes(~drop);

f = fcn_rickplot([2,2,5,5]);
scatter(coor_tsne(:,1),coor_tsne(:,2),fcn_sz(num_nodes,[5,100]),cols_cent(clu_relabel,:),'filled')
fcn_axlims(coor_tsne(:,1),coor_tsne(:,2))

%% Figure 7d

cmap_center = [0 0 1; 1 1 1; 1 0 0];
cmap_center = interp1([0,1,2],cmap_center,linspace(0,2,256));
load ../data/coor_lims.mat
meshdata = load('../data/mesh_data_neuropil.mat');

wid = 6;
x = mx - mn;
x = x/max(x);
dims = [2,2,wid,wid*x(2)];

tgtnames = {'internal_density','local_efficiency','arbour_volume','inputoutput_synapse_similarity','jsd_input2output_synapse_weighted','diversity'};
for t = 1:length(tgtnames)

    tdx = contains(N,tgtnames{t}) & cellfun(@length,N) == length(tgtnames{t});
    vals = z(~drop,tdx);

    prct = prctile(vals,90);
    vec = z(ci(:,1),tdx);


    f = fcn_rickplot(dims);
    ax = axes; hold(ax,'on');
    th = zeros(1,length(meshdata.f));
    for k = 1:length(meshdata.f)
        th(k) = trisurf(meshdata.f{k},meshdata.v{k}(:,1),-meshdata.v{k}(:,2),meshdata.v{k}(:,3));
    end
    axis image off;
    set(th,'edgecolor','none','facecolor',ones(1,3)*0.5,'facealpha',0.05);
    material dull;
    camlight headlight

    scatter3(coor(:,1),-coor(:,2),coor(:,3),5,vec,'filled');
    set(gca,'clim',[-prct,prct]);
    colormap(cmap_center)

    f = fcn_rickplot([2,2,5,5]);
    scatter(coor_tsne(:,1),coor_tsne(:,2),fcn_sz(num_nodes,[5,100]),vals,'filled');
    colormap(cmap_center);
    set(gca,'clim',[-max(abs(vals)),max(abs(vals))]*0.5);

    vec = z(ci(:,1),tdx);

end

%%

z_cent = fcn_dummyvar(clu_relabel,[],1)'*z(~drop,:);
rho_cent = corr(z_cent');

m = dummyvar(clu_relabel)'*dummyvar(clu_relabel);

%% Fig. 8a

c = cidown(~drop,max(cidown) > 3);

z_lodim = z(~drop,:);
rho = (corr(z_lodim'));
vals = [];
for l = 1:size(c,2)
    clu = fcn_dummyvar(c(:,l),[],1);
    r = clu'*(rho*clu);
    vals = [vals; diag(r),ones(length(r),1)*l];
end

nrand = 100;
mj = zeros(nrand,max(vals(:,2)));

for irand = 1:nrand
    disp(irand)
    valsr = [];
    rr = randperm(size(z_lodim,1));
    for l = 1:size(c,2)
        clu = fcn_dummyvar(c(rr,l),[],1);
        r = clu'*(rho*clu);
        valsr = [valsr; diag(r),ones(length(r),1)*l];
    end

    for j = 1:max(valsr(:,2))
        jdx = valsr(:,2) == j;
        mj(irand,j) = nanmedian(valsr(jdx,1));
    end

end

x = vals(:,2) + unifrnd(-0.25,0.25,size(vals(:,2)));
y = vals(:,1);

f = fcn_rickplot([2,2,3,2]);
hold(gca,'on');
plot(1:9,mj,'k');
scatter(x,y,10,'filled');

[xx,yy] = fcn_get_boxplot_skeleton(vals(:,1),vals(:,2));
for i = 1:length(xx)
    plot(xx{i},yy{i})
end
fcn_axlims(x,y)

f = fcn_boxpts(vals(:,1),vals(:,2),fcn_cmap(cmap_center,length(unique(vals(:,2)))));

%% Fig. 8b

meta_clu = fcn_dummyvar(clu_relabel);
vals = [];
for l = 1:size(c,2)

    clu = fcn_dummyvar(c(:,l));
    a = meta_clu'*clu;
    prob = bsxfun(@rdivide,a,sum(a));
    ent = -nansum(prob.*log2(prob),1);

    clu = fcn_dummyvar(c(:,l),[],1);
    aaa = clu'*(rho*clu);

    vals = [vals; ent(:),ones(size(ent(:)))*l,diag(aaa)];
    
    val = ent(c(:,l));
    vec = nan(ncoarse,1);
    vec(~drop) = val;
    vec_neuron = vec(ci(:,1));

    t = vec_neuron;

end

mj = zeros(nrand,9);
for irand = 1:nrand

    disp(irand);

    rr = randperm(size(c,1));
    valsr = [];
    for l = 1:size(c,2)


        clu = fcn_dummyvar(c(rr,l));
        a = meta_clu'*clu;
        prob = bsxfun(@rdivide,a,sum(a));
        ent = -nansum(prob.*log2(prob),1);

        clu = fcn_dummyvar(c(:,l),[],1);
        aaa = clu'*(rho*clu);

        valsr = [valsr; ent(:),ones(size(ent(:)))*l,diag(aaa)];

    end

    for j = 1:max(valsr(:,2))
        jdx = valsr(:,2) == j;
        mj(irand,j) = mean(valsr(jdx,1));
    end

end

x = vals(:,2) + unifrnd(-0.25,0.25,size(vals(:,2)));
y = vals(:,1);

f = fcn_rickplot([2,2,3,2]);
hold(gca,'on');
plot(1:9,mj,'k')
scatter(x,y,10,'filled');

[xx,yy] = fcn_get_boxplot_skeleton(vals(:,1),vals(:,2));
for i = 1:length(xx)
    plot(xx{i},yy{i})
end
fcn_axlims(x,[y; mj(:)])

%% 8c

f = fcn_rickplot([2,2,2,2]);
scatter(vals(:,1),vals(:,3),10,vals(:,2),'o','filled');
fcn_axlims(vals(:,1),vals(:,3))
text(max(vals(:,1)),max(vals(:,3)),sprintf('r=%.2f',corr(vals(:,1),vals(:,3),'rows','pairwise')),...
    'horizontalalignment','right',...
    'verticalalignment','top')

%%

list = cell(1,3);
list{1} = [9,3,19,13,12,8,36,6,7,10,4,33,17,16,2,32];
list{2} = [1,30,25,27,24,34,26,31,22,29,14,28,35,41];
list{3} = [5,21,23,11,37,43];
list{4} = [40,44,45];
list{5} = [42,39,18,20,38,15];
ll = zeros(max(clu_relabel),1);
ll(list{1}) = 1;
ll(list{2}) = 2;

mm = zeros(size(clu_relabel));
for k = 1:2
    for i = 1:length(list{k})
        mm(clu_relabel == list{k}(i)) = k;
    end
end

zz = full(z(~drop,:));
zzz = zz(mm == 1 | mm == 2,:);
nnn = mm(mm == 1 | mm == 2);
pp = zeros(size(zz,2),1);
cohensd = pp;
for i = 1:size(zz,2)
    [~,pp(i)] = ttest2(zz(mm == 1,i),zz(mm == 2,i));
    cohensd(i) = fcn_cohens_d(zzz(:,i),nnn);
end

cohenscols = interp1(linspace(-max(abs(cohensd)),max(abs(cohensd)),256),cmap_center,cohensd);
[~,ii] = sort(cohensd,'descend');
f = fcn_rickplot([2,2,10,6]);
hold(gca,'on');
for j = 1:length(ii)
    bar(j,cohensd(ii(j)),'facecolor',cohenscols(ii(j),:))
end
padj = fcn_linear_step_up(pp,0.01);
sig = find(pp(ii) < padj);
plot(sig,ones(size(sig))*max(abs(cohensd)),'k*')
set(gca,'xtick',1:length(cohensd),'xticklabel',strrep(N(ii),'_','-'))

%% Fig 8f

btotal = zeros(max(clu_relabel));
etotal = btotal;

for l = 1:size(c,2)


    clu = fcn_dummyvar(c(:,l));
    a = meta_clu'*clu;
    b = a*a';

    btotal = btotal + b;

    d = diag(diag(b).^-0.5);
    e = d*(b*d);

    etotal = etotal + e;

end

f = fcn_rickplot([2,2,2,2]);
imagesc(etotal,[-9,9])
colormap(cmap_center)

f = fcn_rickplot([2,2,2,2]);
imagesc(rho_cent,[-1,1])
colormap(cmap_center)

%% Fig 8g

etotalnorm = etotal/etotal(1,1);

h = hist(clu_relabel,1:max(clu_relabel));
m = fcn_mst_plus(etotalnorm,0.4);
f = fcn_rickplot([2,2,5,5]);
g = graph(m);
plot(g,'nodecolor',cols_cent,'markersize',fcn_sz(h,[1,25]))

list = cell(1,3);
list{1} = [9,3,19,13,12,8,36,6,7,10,4,33,17,16,2,32];
list{2} = [1,30,25,27,24,34,26,31,22,29,14,28,35,41];
list{3} = [5,21,23,11,37,43];
list{4} = [40,44,45];
list{5} = [42,39,18,20,38,15];