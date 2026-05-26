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
load ../outputs/community_info.mat
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
%%
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
                    % px = bsxfun(@rdivide,x,sum(x,1));
                    % entx = -nansum(px.*log2(px),1)/log2(length(px));
                    % x = entx;
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
mu = nanmean(Xa);
sd = nanstd(Xa);
z = bsxfun(@minus,Xa,mu);
z = bsxfun(@rdivide,z,sd);
z = full(z);

%%
cc = cidown(orderagreement,:);
commcoor = zeros(size(z,1),4);
count = 0;
for i = 1:max(Xa(:,end))
    for j = 1:max(cc(:,i))
        count = count + 1;
        jdx = find(cc(:,i) == j);
        commcoor(count,1) = i;
        commcoor(count,2) = mean(jdx);
        commcoor(count,3) = j;
        commcoor(count,4) = i;
    end
end
%%
cmap = fcn_cmap([0 0 1; 1 1 1; 1 0 0],256);
R = [];
G = [];
XY = [];
LV = [];
ID = [];
count = 0;
nrand = 100;
clear D

for i = 1:8
    zi = z(Xa(:,end) == i,1:end - 1);
    zp = z(Xa(:,end) == (i + 1),1:end - 1);
    rip = zeros(max(cc(:,i)),1);
    ripeuc = zeros(max(cc(:,i)),1);
    ptrue = rip;
    id = zeros(max(cc(:,i)),2);
    csz = zeros(max(cc(:,i)),3);
    for j = 1:max(cc(:,i))
        count = count + 1;
        jdx = find(cc(:,i) == j);
        yc = mean(jdx);
        p = unique(cc(jdx,i + 1));
        ptrue(j) = p;
        yp = mean(find(cc(:,i + 1) == p));
        XY = [XY; i,yc; i + 1,yp; nan,nan];
        fc = zi(j,:)';
        fp = zp(p,:)';
        dff = fc - fp;
        keep = ~isnan(fc) & ~isnan(fp);
        rip(j) = fcn_lins_concordance(fc(keep),fp(keep));
        ripeuc(j) = nansum(dff.^2).^0.5;
        LV = [LV; count*ones(2,1); nan];
        id(j,:) = [j,p];
        csz(j,1:2) = [sum(ci(:,i) == j),sum(ci(:,i + 1) == p)];
        csz(j,1:3) = sum((ci(:,i) == j) & (ci(:,i + 1) == p))/sum((ci(:,i) == j) | (ci(:,i + 1) == p));
    end
    
    rho = corr(zi',zp','rows','pairwise');
    [~,ppred] = sort(rho','descend');
    match = cumsum(bsxfun(@eq,ppred,ptrue'),1);
    R = [R; rip,ones(size(rip))*i];
    ID = [ID; id];
    M = ripeuc;
    C = M(ci(:,i))';
    D(:,i) = C;
    sz = hist(ci(:,i),1:max(ci(:,i)));
    wt = sz/sum(sz);
    G = [G; ripeuc,ones(size(ripeuc))*i,sz(:),csz];

end

%% Fig 8d
linewid = [0.01,1];
lim = 1;
Rclip = fcn_clip(R(:,1),[0,lim]);
col = interp1([0,lim/2,lim],[0 0 1; 1 1 1; 1 0 0],Rclip);
f = fcn_rickplot([2,2,2,12]);
hold(gca,'on')
wid = interp1([0,lim],linewid,Rclip);
for i = 1:count
    idx = LV == i;
    x = XY(idx,1);
    y = XY(idx,2);
    plot(x,y,'color',col(i,:),'linewidth',wid(i));
end

avg = fcn_dummyvar(G(:,2),[],1)'*G(:,1);
avgcol = interp1([0,max(avg)/2,max(avg)],[0 0 1; 1 1 1; 1 0 0],avg);
f = fcn_rickplot([2,2,2,2]);
hold(gca,'on');
for j = 1:length(avg)
    bar(j,avg(j),'facecolor',avgcol(j,:));
end
%% Fig 8e
clear sd
d = zeros(max(Xa(:,end)),size(z,2) - 1);
for i = 1:max(Xa(:,end))
    idx = Xa(:,end) == i;
    vals = z(idx,1:end - 1);
    sd(i,:) = range(vals);
end
zsd = zscore(sd);
rr = corr(zsd,(1:9)');
jj = find(rr > 0.8);
nn = N(jj);
f = fcn_rickplot([2,2,2,2]);
ph = plot(1:9,zsd,1:9,nanmean(zsd,2));
set(ph(end),'color','k','linewidth',2);
col = interp1(linspace(-1,1,256),cmap,rr);
for i = 1:length(rr)
    set(ph(i),'color',col(i,:))
end
zsdend = zsd(end,jj);
text(9.1*ones(size(jj)),zsdend,strrep(nn,'_','-'))
fcn_axlims(1:9,zsd(:));