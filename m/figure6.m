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
addpath fcns

%%

% categories of community features
%
%   1. connectivity
%   2. space
%   3. annotational
%   4. neurotransmitter

allnodes = 1:n;

a = zeros(1000000,1); acount = 0;
d = a; dcount = 0;
c = a; ccount = 0;
p = a; pcount = 0;

degr = sum(W.TOT,2);
degc = sum(W.TOT,1);
twom = sum(degr);

degrbin = sum(W.TOT > 0,2);
degcbin = sum(W.TOT > 0,1);


annotnames = {'flow','super_class','class','sub_class','cell_type'};
ntnames = fieldnames(W); ntnames(contains(ntnames,'TOT')) = [];

skeletondir = '/Users/rbetzel/Brain networks lab Dropbox/bnbl_main/data/drosophila+reassembled/skeleton_matrices_princeton/';

Xa = [];
for lvl = 1:9
    filename = sprintf('../outputs/community_features.lvl=%i.mat',lvl);
    check = dir(filename);
    if isempty(check)

        ncoarse = max(ci(:,lvl));
        assort = zeros(ncoarse);
        disassort = zeros(ncoarse);
        coreperiph = zeros(ncoarse);
        emptymotif = zeros(ncoarse);
        ties = zeros(ncoarse);

        features = struct();
        for s = 1:max(ci(:,lvl))

            %% update
            fprintf('lvl %i, community %i/%i\n',lvl,s,max(ci(:,lvl)))

            %% 0. get indices of nodes in community and community subgraph
            idx = find(ci(:,lvl) == s);
            w = W.TOT(idx,idx);
            nw = length(w);

            %% 1. connectivity-based features

            % modularity contribution
            pmat = degr(idx)*degc(idx)/twom;
            b = (w - pmat)/twom;
            features.connectivity(s).modularity_contribution = sum(b(:));
            features.connectivity(s).modularity_contribution_normalized = sum(b(:))/(nw^2);

            features.connectivity(s).mean_degree_row = mean(degrbin(idx));
            features.connectivity(s).mean_degree_col = mean(degcbin(idx));
            features.connectivity(s).mean_degree_combined = mean(degcbin(idx) + degrbin(idx)');

            features.connectivity(s).variance_degree_row = std(degrbin(idx));
            features.connectivity(s).variance_degree_col = std(degcbin(idx));
            features.connectivity(s).variance_degree_combined = std(degcbin(idx) + degrbin(idx)');

            features.connectivity(s).skew_degree_row = skewness(degrbin(idx));
            features.connectivity(s).skew_degree_col = skewness(degcbin(idx));
            features.connectivity(s).skew_degree_combined = skewness(degcbin(idx) + degrbin(idx)');

            features.connectivity(s).mean_strength_row = mean(degr(idx));
            features.connectivity(s).mean_strength_col = mean(degc(idx));
            features.connectivity(s).mean_strength_combined = mean(degc(idx) + degr(idx)');

            features.connectivity(s).variance_strength_row = std(degr(idx));
            features.connectivity(s).variance_strength_col = std(degc(idx));
            features.connectivity(s).variance_strength_combined = std(degc(idx) + degr(idx)');

            features.connectivity(s).skew_strength_row = skewness(degr(idx));
            features.connectivity(s).skew_strength_col = skewness(degc(idx));
            features.connectivity(s).skew_strength_combined = skewness(degc(idx) + degr(idx)');

            % size = number of nodes
            features.connectivity(s).size = nw;

            % edges = number of internal edges
            features.connectivity(s).edges = nnz(w);

            % volume = internal number of synapses
            features.connectivity(s).volume = sum(sum(w,2) + sum(w,1)');

            % strength_asymmetry = RMS difference between in/out strength
            features.connectivity(s).strength_asymmetry = mean((sum(w,2) - sum(w,1)').^2).^0.5;

            % degree_asymmetry = RMS difference between in/out degree
            features.connectivity(s).strength_asymmetry = mean((sum(w > 0,2) - sum(w > 0,1)').^2).^0.5; % need to update this, but let's not change it now

            % internal_density = normalized number of edges
            features.connectivity(s).internal_density = features.connectivity(s).edges/(nw*(nw - 1));

            % normalized_volume = normalized volume
            features.connectivity(s).normalized_volume = features.connectivity(s).volume/(nw*(nw - 1));

            % local efficiency
            g = digraph(w > 0);
            dist = distances(g);
            e = 1./dist;
            e(1:(nw + 1):end) = nan;

            % local efficiency = mean 1/path_length
            features.connectivity(s).local_efficiency = nanmean(e(:));

            % local efficiency = max path length (for pairs that are reachable)
            features.connectivity(s).diameter = max(dist(~isinf(dist)));

            dist(1:(nw + 1):end) = nan;

            % mean_reachability = mean fraction of nodes that are reachable from
            % any other node
            features.connectivity(s).mean_reachability = mean(nanmean(~isinf(dist),2));

            jdx = setdiff(allnodes,idx);

            vs  = sum(sum(W.TOT(idx,:))) + sum(sum(W.TOT(:,idx)));
            vsc = sum(sum(W.TOT(jdx,:))) + sum(sum(W.TOT(:,jdx)));

            cut_weight = sum(sum(W.TOT(idx,jdx))) + sum(sum(W.TOT(jdx,idx)));

            % cut_weight = total weight from community to other nodes
            features.connectivity(s).cut_weight = cut_weight;

            % normalized_cut_weight = total weight from community to other nodes
            % normalized by number of possible edges
            features.connectivity(s).normalized_cut_weight = cut_weight/(2*length(idx)*length(jdx));

            % conductance = total weight from community to other nodes
            features.connectivity(s).conductance = cut_weight/min(vs,vsc);

            % expansion = total weight from community to other nodes
            features.connectivity(s).expansion = cut_weight/length(idx);

            % boundary_cut =
            int_weight = sum(sum(w));
            features.connectivity(s).boundary_fraction = cut_weight/(cut_weight + int_weight);

            % input2output_ratio = ratio of incoming number of synapses to number
            % of outgoing synapses
            wi = W.TOT(:,idx); wi(idx,:) = 0; si = sum(wi,2); ki = sum(wi > 0,2);
            wo = W.TOT(idx,:); wo(:,idx) = 0; so = sum(wo,1); ko = sum(wo > 0,1);
            features.connectivity(s).input2output_synapse_ratio = full(sum(si)/sum(so));
            features.connectivity(s).input2output_edge_ratio = full(sum(ki)/sum(ko));

            features.connectivity(s).inputoutput_synapse_similarity = 1 - pdist2(full(si'),full(so),'cosine');
            features.connectivity(s).inputoutput_edge_similarity = 1 - pdist2(full(ki'),full(ko),'cosine');
            features.connectivity(s).inputoutput_overlap = length(intersect(find(full(ki)),find(full(ko))))/length(union(find(full(ki)),find(full(ko))));

            vals = 1:ncoarse;
            vals(s) = [];
            features.connectivity(s).segregation_max = synapse_count_norm(s,s) - max(max(synapse_count_norm(s,vals)),max(synapse_count_norm(vals,s)));
            features.connectivity(s).segregation_mean = synapse_count_norm(s,s) - (mean(synapse_count_norm(s,vals)) + mean(synapse_count_norm(vals,s)))/2;

            % assortativity measures (assort, core, periphery, disassort)
            for t = 1:ncoarse
                if s ~= t
                    aa = synapse_count_norm(s,s);
                    bb = synapse_count_norm(t,t);
                    ab = synapse_count_norm(s,t);
                    if min(aa,bb) > ab
                        assort(s,t) = 1;
                        acount = acount + 1;
                        a(acount) = s;
                        acount = acount + 1;
                        a(acount) = t;
                    elseif ab > max(aa,bb)
                        disassort(s,t) = 1;
                        dcount = dcount + 1;
                        d(dcount) = s;
                        dcount = dcount + 1;
                        d(dcount) = t;
                    elseif (aa > ab) & (ab > bb)
                        coreperiph(s,t) = 1;
                        ccount = ccount + 1;
                        c(ccount) = s;
                        pcount = pcount + 1;
                        p(pcount) = t;
                    elseif (bb > ab) & (ab > aa)
                        coreperiph(s,t) = 1;
                        ccount = ccount + 1;
                        c(ccount) = t;
                        pcount = pcount + 1;
                        p(pcount) = s;
                    elseif sum([aa,bb,ab] == 0) > 1
                        emptymotif(s,t) = 1;
                    elseif length(unique([aa,bb,ab])) < 3
                        ties(s,t) = 1;
                    end
                end
            end

            %% 2. space-based features

            ids = root_id(idx);
            skeletonxyz = zeros(1000000,3);
            scount = 0;
            total_vol = zeros(length(ids),1);
            total_len = total_vol;
            num_branches = total_len;
            for j = 1:length(ids)
                tmp = load(fullfile(skeletondir,sprintf('%s.mat',ids{j})));
                jdx = (scount + 1):(scount + size(tmp.coor,1));
                scount = scount + size(tmp.coor,1);
                skeletonxyz(jdx,:) = tmp.coor;
                [stubu,stubv,thickness] = find(tmp.a);
                len = sum((tmp.coor(stubu,:) - tmp.coor(stubv,:)).^2,2).^0.5;
                vol = len(:).*thickness(:);
                total_vol(j) = sum(vol);
                total_len(j) = sum(len);
                num_branches(j) = sum(sum(tmp.a > 0,2) > 1);
            end
            skeletonxyz = skeletonxyz(1:scount,:);

            xyz = coor(idx,:);
            soma_centroid = mean(xyz,1);
            arbour_centroid = nanmean(skeletonxyz,1);
            features.spatial(s).soma_centroid = soma_centroid';
            features.spatial(s).arbour_centroid = arbour_centroid';

            diffs = skeletonxyz - soma_centroid;
            dists = sqrt(sum(diffs.^2,2));
            features.spatial(s).radius_of_gyration = sqrt(mean(dists.^2));

            pd = pdist(xyz);
            features.spatial(s).mean_pairwise_soma2soma_distance = mean(pd);

            [verts, hull_vol] = convhull(skeletonxyz(:,1),skeletonxyz(:,2),skeletonxyz(:,3));
            features.spatial(s).convex_hull_vol = hull_vol;
            sphere_volume = (4/3)*pi*features.spatial(s).radius_of_gyration^3;
            features.spatial(s).compactness_ratio = hull_vol/sphere_volume;
            features.spatial(s).mean_num_branches = mean(num_branches);
            features.spatial(s).arbour_length = sum(total_len);
            features.spatial(s).arbour_volume = sum(total_vol);

            %% 3. annotational features

            si = find(sum(wi,2) > 0);
            so = find(sum(wo,1) > 0);
            sr = idx;

            swi = sum(wi,2); swi = swi(si);
            swo = sum(wo,1); swo = swo(so);
            swr = sum(w,1) + sum(w,2)';

            jsd_io_wei = zeros(size(annotnames));
            jsd_ir_wei = zeros(size(annotnames));
            jsd_ro_wei = zeros(size(annotnames));
            jsd_io = zeros(size(annotnames));
            jsd_ir = zeros(size(annotnames));
            jsd_ro = zeros(size(annotnames));
            dff = jsd_ro;
            dff_wei = dff;
            for j = 1:length(annotnames)

                ncats = max(labels.(annotnames{j}));
                valsi = labels.(annotnames{j})(si); hi = hist(valsi,1:ncats)/length(valsi);
                valso = labels.(annotnames{j})(so); ho = hist(valso,1:ncats)/length(valso);
                valsr = labels.(annotnames{j})(sr); hr = hist(valsr,1:ncats)/length(valsr);

                hiw = hi;
                how = ho;
                hrw = hr;
                for k = 1:ncats

                    kdxi = valsi == k;
                    kdxo = valso == k;
                    kdxr = valsr == k;

                    hiw(k) = sum(swi(kdxi));
                    how(k) = sum(swo(kdxo));
                    hrw(k) = sum(swr(kdxr));

                end
                hiw = hiw/sum(hiw);
                how = how/sum(how);
                hrw = hrw/sum(hrw);

                jsd_io(j) = 1 - fcn_jsd(hi,ho);
                jsd_ir(j) = 1 - fcn_jsd(hi,hr);
                jsd_ro(j) = 1 - fcn_jsd(hr,ho);

                enti = -nansum(hi.*log2(hi))/log2(ncats);
                ento = -nansum(ho.*log2(ho))/log2(ncats);

                dff(j) = enti - ento;

                jsd_io_wei(j) = 1 - fcn_jsd(hiw,how);
                jsd_ir_wei(j) = 1 - fcn_jsd(hiw,hrw);
                jsd_ro_wei(j) = 1 - fcn_jsd(hrw,how);

                enti_wei = -nansum(hiw.*log2(hiw))/log2(ncats);
                ento_wei = -nansum(how.*log2(how))/log2(ncats);

                dff_wei(j) = enti_wei - ento_wei;

            end

            features.annotations(s).jsd_input2output = mean(jsd_io);
            features.annotations(s).jsd_input2community = mean(jsd_ir);
            features.annotations(s).jsd_community2output = mean(jsd_ro);
            features.annotations(s).delta_annotation_entropy = mean(dff);

            features.annotations(s).jsd_input2output_synapse_weighted = mean(jsd_io_wei);
            features.annotations(s).jsd_input2community_synapse_weighted = mean(jsd_ir_wei);
            features.annotations(s).jsd_community2output_synapse_weighted = mean(jsd_ro_wei);
            features.annotations(s).delta_annotation_entropy_synapse_weighted = mean(dff_wei);

            %% 4. neurotransmitter features

            jdx = setdiff(allnodes,idx);

            for nt = 1:length(ntnames)

                Wnt = full(W.(ntnames{nt})(idx,idx));

                features.neurotransmitter(s).within_community_edge_count(nt) = nnz(Wnt);
                features.neurotransmitter(s).within_community_synapse_count(nt) = sum(Wnt(:));

                features.neurotransmitter(s).within_community_edge_density(nt) = features.neurotransmitter(s).within_community_edge_count(nt)/(length(idx)*(length(idx) - 1));
                features.neurotransmitter(s).within_community_synapse_density(nt) = features.neurotransmitter(s).within_community_synapse_count(nt)/(length(idx)*(length(idx) - 1));

                Wnt = full(W.(ntnames{nt})(idx,jdx));

                features.neurotransmitter(s).outgoing_community_edge_count(nt) = nnz(Wnt);
                features.neurotransmitter(s).outgoing_community_synapse_count(nt) = sum(Wnt(:));

                features.neurotransmitter(s).outgoing_community_edge_density(nt) = features.neurotransmitter(s).within_community_edge_count(nt)/(length(idx)*(length(idx) - 1));
                features.neurotransmitter(s).outgoing_community_synapse_density(nt) = features.neurotransmitter(s).within_community_synapse_count(nt)/(length(idx)*(length(idx) - 1));

                Wnt = full(W.(ntnames{nt})(jdx,idx));

                features.neurotransmitter(s).incoming_community_edge_count(nt) = nnz(Wnt);
                features.neurotransmitter(s).incoming_community_synapse_count(nt) = sum(Wnt(:));

                features.neurotransmitter(s).incoming_community_edge_density(nt) = features.neurotransmitter(s).within_community_edge_count(nt)/(length(idx)*(length(idx) - 1));
                features.neurotransmitter(s).incoming_community_synapse_density(nt) = features.neurotransmitter(s).within_community_synapse_count(nt)/(length(idx)*(length(idx) - 1));

            end

            %% 5. SBM-specific features


            if lvl == 1
                % for k = 1:size(ci,2)
                for k = 2:9
                    parent = unique(ci(idx,k));
                    pdx = find(ci(:,k) == parent);
                    features.sbm(s).persistence(k) = length(intersect(pdx,idx))/length(union(pdx,idx));
                    wk = W.TOT(pdx,pdx);
                    edgesk = nnz(wk);
                    nkw = length(pdx);
                    densityk = edgesk/(nkw*(nkw - 1));
                    features.sbm(s).difference_in_density(k) = features.connectivity(s).internal_density - densityk;
                end
                features.sbm(s).persistence(1) = nan;
                features.sbm(s).difference_in_density(1) = nan;
            end


        end
        %%
        a = a(1:acount); d = d(1:dcount); c = c(1:ccount); p = p(1:pcount);
        aa = hist(a,1:ncoarse); dd = hist(d,1:ncoarse); cc = hist(c,1:ncoarse); pp = hist(p,1:ncoarse);
        hh = [aa',dd',cc',pp'];
        qq = bsxfun(@rdivide,hh,sum(hh,2));
        ee = -nansum(qq.*log2(qq),2);

        assort_data = [hh,ee];

        save(filename,'features','assort_data')
    else
        load(filename);
        %%

        ncoarse = max(ci(:,lvl));

        features.spatial = rmfield(features.spatial,'soma_centroid');
        features.spatial = rmfield(features.spatial,'arbour_centroid');

        fn = fieldnames(features);
        fn(contains(fn,'neurotransmitter')) = [];
        fn(contains(fn,'sbm')) = [];

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
                end
                if contains(fn{j},'sbm')
                    m = length(x)/ncoarse;
                    x = reshape(x,[m,ncoarse])';
                    x = nanmean(x,2);
                end
                [nr,nc] = size(x);
                if nc == ncoarse
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


        N = [N; {'assortativity'}; {'disassortativity'}; {'core'}; {'periphery'}; {'diversity'}];

        Xa = [Xa; X, ones(size(X,1),1)*lvl];

        %%
    end
end
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

Xr = X;
z = bsxfun(@rdivide,bsxfun(@minus,Xr,nanmean(Xr)),nanstd(Xr));
N = [N; {'assortativity'}; {'disassortativity'}; {'core'}; {'periphery'}; {'diversity'}];

drop = isnan(sum(full(z),2));

[coeff,score,latent] = pca(full(z));
latent = latent/sum(latent);
clatent = cumsum(latent);
varexp = 0.85;
feat_lodim = score(~drop,1:find(clatent >= varexp,1,'first'));

coeff_lodim = coeff(:,1:size(feat_lodim,2));
[~,ii] = max(abs(coeff_lodim),[],2);
jj = [];
for i = 1:size(coeff_lodim,2)
    kk = find(ii == i);
    vv = coeff_lodim(kk,i);
    ss = sign(vv);
    jj = [jj; kk(ss > 0); kk(ss < 0)];
end

cmap_center = interp1([0,1,2],[0 0 1; 1 1 1; 1 0 0],linspace(0,2,256));

%% Fig 6b
f = fcn_rickplot([2,2,5,5]);
imagesc(coeff_lodim(jj,:),[-1,1]*0.5)
colormap(cmap_center);
set(gca,'ytick',1:size(z,2),'yticklabel',strrep(N(jj),'_','-'))

%% Fig 6d

markersz = 5;
for i = 1:3
    vec = score(ci(:,1),i);

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

    scatter3(coor(:,1),-coor(:,2),coor(:,3),markersz*(abs(vec) + eps)/max(abs(vec)),vec,'filled');
    set(gca,'clim',[-5,5])
    colormap(cmap_center)
    axis image;
    view([0,90])

end


%% Fig 6f
[~,ii] = max(abs(feat_lodim),[],2);
jj = [];
for i = 1:size(feat_lodim,2)
    kk = find(ii == i);
    vv = feat_lodim(kk,i);
    ss = sign(vv);
    jj = [jj; kk(ss > 0); kk(ss < 0)];
end

scale_factor = 5;
f = fcn_rickplot([2,2,2,4]);
imagesc(feat_lodim(jj,:),[-1,1]*scale_factor)
colormap(cmap_center);