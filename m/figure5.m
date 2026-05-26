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
%% Figure 5a,b,c,d
gts = fieldnames(names);
for t = 1:4

    L = dummyvar(labels.(tgts{t}));
    drop = contains(names.(tgts{t}),' ') | cellfun(@isempty,names.(tgts{t}));
    nm = names.(tgts{t})(~drop);
    L = L(:,~drop);
    Ai = zeros(size(L,2));
    n = size(L,1);
    C = dummyvar(ci(:,1));
    LC = C'*L;
    M = LC'*LC;

    filename = sprintf('../outputs/order_coassignment.%s.mat',tgts{t});
    check = dir(filename);
    if isempty(check)
        idxsort = 1:length(M);
        bb = inf;
        for irep = 1:100
            fprintf('%i,%i\n',t,irep)
            [~,ii,cc] = reorder_matrix(M,'line',1,100000);
            if cc < bb
                bb = cc;
                idxsort = ii;
            end
        end

        save(filename,'M','idxsort')

    else

        load(filename);
    end

    [u,v,w] = find(triu(M(idxsort,idxsort)));
    edges = (v - 1)*length(M) + u;

    mn = min(log10(w));
    mx = max(log10(w));
    cols = interp1(linspace(mn,mx,size(cmap,1)),cmap,log10(w));

    g = get_components(M);

    lb = labels.(tgts{t});
    if sum(drop)
    lb_drop = lb; lb_drop(lb_drop == find(drop)) = 0;
    gl = zeros(size(lb_drop));
    gl(lb_drop > 0) = g(lb_drop(lb_drop > 0) - sum(drop));
    else
        gl = g(lb);
    end

    markersz = 2.5;
    f = fcn_rickplot([2,2,5,5]);
    imagesc(log10(M(idxsort,idxsort)))
    set(gca,...
        'xtick',1:length(idxsort),'xticklabel',strrep(nm(idxsort),'_','-'),...
        'ytick',1:length(idxsort),'yticklabel',strrep(nm(idxsort),'_','-'));
    
    [u,v,w] = find(M(idxsort,idxsort));
    edges = (v - 1)*length(M) + u;

    mn = min(log10(w));
    mx = max(log10(w));
    cmap = parula(256);
    cols = interp1(linspace(mn,mx,size(cmap,1)),cmap,log10(w));

    %%

    g = graph(M(idxsort,idxsort));
    figure; ph = plot(g);
    xy = [ph.XData',ph.YData'];
    close(gcf);
    
    
    deg = sum(M(idxsort,idxsort));


    wid = fcn_sz(log10(w),[0.5,3]);
    f = fcn_rickplot([2,2,5,5]);
    ax = axes; hold(ax,'on');
    for j = 1:length(u)
        x = [xy(u(j),1),xy(v(j),1)];
        y = [xy(u(j),2),xy(v(j),2)];
        plot(x,y,'linewidth',wid(j),'color',cols(j,:));
        
    end
    scatter(xy(:,1),xy(:,2),fcn_sz(log10(diag(M(idxsort,idxsort))),[15,200]),interp1(linspace(mn,mx,size(cmap,1)),cmap,log10(diag(M(idxsort,idxsort)))),'filled')
    fcn_axlims(xy(:,1),xy(:,2));
    text(xy(:,1),xy(:,2),strrep(nm(idxsort),'_','-'));
    axis image off;
    

end