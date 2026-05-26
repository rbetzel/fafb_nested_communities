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

%%

% calculate community motif counts for all hierarchical levels
mat = cell(1,9);
for ilevel = 1:9
    sc = W.TOT;
    c = dummyvar(ci(:,ilevel));
    den = sum(c)'*sum(c);
    synapse_count = c'*(sc*c);
    synapse_count_norm = synapse_count./den;

    ncoarse = length(synapse_count_norm);
    assort = nan(ncoarse);
    disassort = nan(ncoarse);
    coreperiph = nan(ncoarse);
    emptymotif = nan(ncoarse);
    ties = nan(ncoarse);
    a = zeros(1000000,1); acount = 0;
    d = a; dcount = 0;
    c = a; ccount = 0;
    p = a; pcount = 0;
    ecount = 0;
    tcount = 0;
    %%
    for i = 1:length(synapse_count_norm)
        for j = 1:length(synapse_count_norm)
            if i ~= j
                aa = synapse_count_norm(i,i);
                bb = synapse_count_norm(j,j);
                ab = synapse_count_norm(i,j);

                % if both within community blocks have greater density than
                % the between-community block, then the motif is
                % assortative
                if min(aa,bb) > ab
                    assort(i,j) = 1;
                    acount = acount + 1;
                    a(acount) = i;
                    acount = acount + 1;
                    a(acount) = j;

                % if the between is greater than both within, then the
                % motif is disassortative
                elseif ab > max(aa,bb)
                    disassort(i,j) = 1;
                    dcount = dcount + 1;
                    d(dcount) = i;
                    dcount = dcount + 1;
                    d(dcount) = j;

                % if a-a density is greater than a-b, and a-b is greater
                % than b-b, then core-periperhy (a the core, b is the
                % periphery)
                elseif (aa > ab) & (ab > bb)
                    coreperiph(i,j) = 1;
                    ccount = ccount + 1;
                    c(ccount) = i;
                    pcount = pcount + 1;
                    p(pcount) = j;

                % if b-b density is greater than a-b, and a-b is greater
                % than a,a, then core-periperhy (b the core, a is the
                % periphery)
                elseif (bb > ab) & (ab > aa)
                    coreperiph(i,j) = 1;
                    ccount = ccount + 1;
                    c(ccount) = j;
                    pcount = pcount + 1;
                    p(pcount) = i;
                
                % if there are multiple zeros, then it's hard to classify
                % motify type
                elseif sum([aa,bb,ab] == 0) > 1
                    emptymotif(i,j) = 1;
                    ecount = ecount + 1;

                % if there are ties, then classification is also hard -
                % these tend to involve at least one small community 
                % (1 or 2 neuron)
                elseif length(unique([aa,bb,ab])) < 3
                    [aa,bb,ab,sum(ci(:,ilevel) == j),sum(ci(:,ilevel) == i)]
                    ties(i,j) = 1;
                    tcount = tcount + 1;
                end
            end
        end
    end

    if ilevel == 1

        %% Figure 3a

        % assortativity matrix
        f = fcn_rickplot([2,2,2,2]);
        imagesc(assort(order,order) == 1);
        colormap([1 1 1; 0 0 0]);

        %% Figure 3b

        % disassortativity matrix
        f = fcn_rickplot([2,2,2,2]);
        imagesc(disassort(order,order) == 1);
        colormap([1 1 1; 0 0 0]);

        %% Figure 3c

        % core-periphery matrix
        f = fcn_rickplot([2,2,2,2]);
        imagesc(coreperiph(order,order) == 1);
        colormap([1 1 1; 0 0 0]);

    end

    % return codes for community interactions:
    % 1 = assort
    % 2 = disassort
    % 3 = core-periphery
    % 4 = empty
    % 5 = tie

    mat{ilevel} = ~isnan(assort) + 2*~isnan(disassort) + 3*~isnan(coreperiph) + 4*~isnan(emptymotif) + 5*~isnan(ties);

end

%%

% calculate relative probability of each motif type at every hierarchical
% level

A = zeros(size(cidown,1),9,5);

cmap = [[46 49 146]/255; [249,237,50]/255; 1 0 0];

c = A;
for ii = 1:9

    total = mat{ii};
    sz = sum(dummyvar(ci(:,ii)));
    den = sz'*sz;
    
    vals = zeros(length(total),5);
    for i = 1:5
        vals(:,i) = sum((total == i).*den,2) + sum((total == i).*den,1)';
    end
    valsprob = bsxfun(@rdivide,vals,sum(vals,2));

    for k = 1:5
        c(:,ii,k) = valsprob(cidown(:,ii),k);
    end

end

%%

% draw motif probabilities for each neuron at each level

cmapnicebase = [[46 49 146]/255; [249,237,50]/255; 1 0 0];
ccc = cmapnicebase([3,1,2],:);
lims = [0.9,1; 0,0.15; 0,0.15];
motiftype = {'assort','disassort','coreperiph'};
for i = 1:3
    cmap = interp1([0,1],[1 1 1; ccc(i,:)],linspace(0,1,256));
    
    %% Figures 3d, 3e, 3f
    f = fcn_rickplot([2,2,3,2]);
    imagesc(c(order,:,i),lims(i,:));
    colormap(cmap);

end

%% Fig. 3g

% draw assortative probabilities in anatomical space

cmap = [0*ones(1,3); [46 49 146]/255; [249,237,50]/255; 1 0 0];
cmap = interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,256));

meshdata = load('../data/mesh_data_neuropil.mat');
load ../data/coor_lims.mat
load ../data/coordinates.mat

x = mx - mn;
x = x/max(x);
wid = 6;
dims = [2,2,wid,wid*x(2)];

total = mat{1};

vals = zeros(length(total),5);
for i = 1:5
    vals(:,i) = sum(total == i,2) + sum(total == i,1)';
end
valsprob = bsxfun(@rdivide,vals,sum(vals,2));

for i = 1

    vals_neurons = valsprob(ci(:,1),i);
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

    scatter3(coor(:,1),-coor(:,2),coor(:,3),fcn_sz(vals_neurons,[0.1,2.5]),(vals_neurons),'filled');
    view([0,90])
    set(gca,'clim',prctile((vals_neurons),[1,99]));
    prctile((vals_neurons),[1,99])
    colormap(cmap);
    axis image off;
 
end

%% Fig. 3h

% draw nonassortative (disassort + core + periph) probabilities in anatomical space

vals_nonassort = sum(valsprob(:,2:3),2);
vals_neurons = vals_nonassort(ci(:,1));

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

scatter3(coor(:,1),-coor(:,2),coor(:,3),fcn_sz(vals_neurons,[0.1,10]),(vals_neurons),'filled');
view([0,90])
set(gca,'clim',prctile((vals_neurons),[1,99]));
prctile((vals_neurons),[1,99])
colormap(cmap);
axis image off;

%% Fig. 3i

% draw entropy over all motif classes in anatomical space

entropy = -nansum(valsprob.*log2(valsprob),2);
vals_neurons = entropy(ci(:,1));

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

scatter3(coor(:,1),-coor(:,2),coor(:,3),fcn_sz(vals_neurons,[0.1,10]),(vals_neurons),'filled');
view([0,90])
set(gca,'clim',prctile((vals_neurons),[1,99]));
prctile((vals_neurons),[1,99])
colormap(cmap);
axis image off;

%% 

% calculate core/periphery counts separately (at finest hierarchical level)

c = dummyvar(ci(:,1));
den = sum(c)'*sum(c);
synapse_count = c'*(sc*c);
synapse_count_norm = synapse_count./den;

nonassort = mat{1} == 3;
[u,v] = find(nonassort);
corecount = zeros(length(nonassort),1);
periphcount = corecount;
for i = 1:length(u)
    aa = synapse_count_norm(u(i),u(i));
    bb = synapse_count_norm(v(i),v(i));
    ab = synapse_count_norm(u(i),v(i));
    if aa > bb
        corecount(u(i)) = corecount(u(i)) + 1;
        periphcount(v(i)) = periphcount(v(i)) + 1;
    else
        corecount(v(i)) = corecount(v(i)) + 1;
        periphcount(u(i)) = periphcount(u(i)) + 1;
    end
end
cfull = corecount(ci(:,1));
pfull = periphcount(ci(:,1));

%% Fig. 3j

% draw core neurons in anatomical space

vals_neurons = (cfull);
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

scatter3(coor(:,1),-coor(:,2),coor(:,3),fcn_sz(vals_neurons,[0.1,10]),vals_neurons,'filled');
view([0,90])
set(gca,'clim',[0,150]);
colormap(cmap);
axis image off;

%% Fig. 3k

% draw perihperal neurons in anatomical space

vals_neurons = (pfull);
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

scatter3(coor(:,1),-coor(:,2),coor(:,3),fcn_sz(vals_neurons,[0.1,10]),vals_neurons,'filled');
view([0,90])
set(gca,'clim',[0,150]);
colormap(cmap);
axis image off;

%% Fig. 3l

% draw most extreme core/periphery neurons in anatomical space

thr = 75;
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

aa = cfull > thr;
bb = pfull > thr;
cc = zeros(size(aa));
cc(aa & ~bb) = 2;
cc(aa & bb) = 0;
cc(~aa & bb) = 1;
mm = aa | bb;

scatter3(coor(:,1),-coor(:,2),coor(:,3),(cc > 0)*2.5 + 0.2,cc,'filled');
view([0,90])
set(gca,'clim',[0,2]);
colormap([ones(1,3)*0.5; [249,237,50]/255; 1 0 0]);

%% Fig. 3m

% draw the different as bar plots

load ../data/annotations.mat

dff = cfull - pfull;
ccc = [0.9765    0.9294    0.1961; 0.75*ones(1,3); 1 0 0];

tgtlist = {'flow','super_class','class','sub_class'};
for t = 1:4

    tgt = tgtlist{t};
    nn = names.(tgt);
    mm = labels.(tgt);

    drop = contains(nn,' ');
    nn(drop) = [];

    ll = dummyvar(labels.(tgt));
    ll = bsxfun(@rdivide,ll,sum(ll,1));

    ll = ll(:,~drop);
    dd = dff'*ll;
    [~,idxsort] = sort(dd,'descend');
    cols = interp1(linspace(-max(abs(dd)),max(abs(dd)),3),ccc,dd);
    
    f = fcn_rickplot([2,2,6,5]);
    ax = axes; hold(ax,'on');
    for j = 1:length(dd)
        bar(j,abs(dd(idxsort(j))),'facecolor',cols(idxsort(j),:));
    end
    set(gca,...
        'xlim',[0,length(dd) + 1],...
        'xtick',1:length(dd),...
        'xticklabel',strrep(nn(idxsort),'_','-'))

end