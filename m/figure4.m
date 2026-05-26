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
tgts = fieldnames(names);
tgts = tgts(1:5);
nrand = 1000;
R = zeros(size(ci,1),nrand,'single');
for i = 1:nrand
    R(:,i) = randperm(size(R,1));
end

%%
for i = 1:9
    annotname = sprintf('../outputs/zannot.lvl=%1.2i.mat',i);
    check = dir(annotname);
    if isempty(check)
        C = dummyvar(ci(:,i));
        for t = 1:length(tgts)
            L = dummyvar(labels.(tgts{t}));
            Or = zeros(size(C,2),size(L,2),nrand);
            O = C'*L;
            for irand = 1:nrand
                Cr = C(R(:,irand),:);
                Or(:,:,irand) = Cr'*L;
                fprintf('hierarchical level = %i, tgt = %s, rand = %i/%i\n',i,tgts{t},irand,nrand);
            end
            o.(tgts{t}) = O;
            m.(tgts{t}) = nanmean(Or,3);
            s.(tgts{t}) = nanstd(Or,[],3);
            z.(tgts{t}) = (O - m.(tgts{t}))./s.(tgts{t});
        end
        save(annotname,'m','s','z','o')
    end
end

%% Fig 4a, 4b

i = 1;
annotname = sprintf('../outputs/zannot.lvl=%1.2i.mat',i);
load(annotname);

fn = fieldnames(o);
for i = 1:4

    nm = names.(fn{i});
    drop = contains(nm,' ') | cellfun(@isempty,nm);
    nm = nm(~drop);

    zmat = z.(fn{i});
    zmat(:,drop) = [];
    zmat(isinf(zmat) | isnan(zmat)) = 0;

    [mx,idxmx] = max(zmat,[],2);

    movetoend = mx <= 1;

    zmat_keep = zmat(~movetoend,:);
    [~,idxsort] = sort(idxmx(~movetoend));

    f = fcn_rickplot([2,2,2,2]);
    imagesc(zmat_keep(idxsort,:),[-10,10]);
    set(gca,'xtick',1:size(zmat_keep,2),'xticklabel',strrep(nm,'_','-'))

    omat = o.(fn{i});
    omat(:,drop) = [];
    omat(isinf(omat) | isnan(omat)) = 0;

    omat_keep = omat(~movetoend,:);

    f = fcn_rickplot([2,2,2,2]);
    imagesc(omat_keep(idxsort,:),[0,max(omat_keep(:))*0.25]);
    set(gca,'xtick',1:size(zmat_keep,2),'xticklabel',strrep(nm,'_','-'));
    colormap(flipud(gray))

end

%% Fig. 4c

ind = dummyvar(ci(:,1));

for i = 1:4

    filename = sprintf('../outputs/dice_overlap.%s.mat',fn{i});
    check = dir(filename);

    nm = names.(fn{i});

    if isempty(check)

        lbl = labels.(fn{i});
        unq = unique(nonzeros(lbl));
        
        M = cell(1,length(unq));
        D = M;

        for l = 1:length(unq)

            close all;

            tgt_annot = +(lbl == unq(l));

            possible = sum(bsxfun(@and,ind,tgt_annot),1) > 0;

            ind_possible = ind(:,possible);

            mask_rep = zeros(10,length(ind_possible));

            dc_rep = zeros(10,1);

            max_count = 1000;
            max_iter = 10000;

            for irep = 1:10

                col = rand(1,3);

                mask = rand(1,size(ind_possible,2)) > 0.5;

                map = sum(ind_possible(:,mask),2);

                dc = dice(map,tgt_annot);

                count = 0;

                iter = 0;

                while (iter < max_iter) & (count < max_count)

                    iter = iter + 1;
                    count = count + 1;

                    r = randi(size(mask));

                    if mask(r) > 0

                        map_iter = map - ind_possible(:,r);

                    else

                        map_iter = map + ind_possible(:,r);

                    end

                    dc_iter = dice(map_iter,tgt_annot);

                    if dc_iter > dc

                        count = 0;

                        dc = dc_iter;

                        mask(r) = ~mask(r);

                        map = map_iter;

                    end

                    if mod(iter,100) == 0
                        plot(iter,dc,'o','color',col);
                        hold on;
                        title(sprintf('%s,%s,%i/%i',fn{i},nm{l},irep,10));
                        drawnow;
                    end

                end

                mask_rep(irep,:) = map;
                dc_rep(irep) = dc;

            end

            M{l} = mask_rep;
            D{l} = dc_rep;

        end

        save(filename,'M','D')

    else

        load(filename);

    end

    %%

    dc_avg = cellfun(@mean,D);

    f = fcn_rickplot([2,2,4,2]);

    bar(dc_avg);
    set(gca,'xtick',1:length(dc_avg),'xticklabel',strrep(nm,'_','-'))

    %%

end

%% Fig 4e

for i = 1:9

    annotname = sprintf('../outputs/zannot.lvl=%1.2i.mat',i);
    load(annotname);
    
    sz = sum(dummyvar(ci(:,i)))';
    sz = sz'/sum(sz);

    fn = fieldnames(z);

    for j = 1:5
        zw = z.(fn{j});
        zw(isinf(zw)) = nan;
        mx = nanmax(zw,[],2);
        zw_avg(i,j) = sum(mx(:).*sz(:));
    end

end

%% figure 4f

zcutoff = 2.2;
V = zeros(length(W.TOT),9);
for i = 1:9
    load(sprintf('../outputs/zannot.lvl=%1.2i.mat',i))
    tgts = fieldnames(names);
    S = zeros(max(ci(:,i)),4);
    for itgt = 1:4
        nm = names.(tgts{itgt});
        drop = contains(nm,' ');
        zz = z.(tgts{itgt});
        zz(isinf(zz) | isnan(zz)) = 0;
        zz = zz(:,~drop);
        s = sum(zz >= zcutoff,2);
        S(:,itgt) = s;
    end
    savg = nanmean(S,2);
    V(:,i) = savg(ci(:,i));
end

meshdata = load('../data/mesh_data_neuropil.mat');
load ../data/coordinates.mat
load ../data/coor_lims.mat

nn = nanmean(zscore(V),2);
idx = nn >= prctile(nn,95);

nlims = 51;
x = mx - mn;
x = x/max(x);
wid = 6;
dims = [2,2,wid,wid*x(2)];
markersz = 2.5;

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
scatter3(coor(~idx,1),-coor(~idx,2),coor(~idx,3),markersz*0.25,nn(~idx,1),'filled')
scatter3(coor(idx,1),-coor(idx,2),coor(idx,3),markersz*5,nn(idx,1),'filled','markeredgecolor','k')
set(gca,'clim',[0,8])