function [e,fitness,fitness_tracker] = fcn_consensus_communities_greedy(c,numiter,max_counter)
% Consensus clustering by greedy ARI maximization with delta updates
%
% c: (n x m) matrix of partitions (n nodes, m ensemble members)
% numiter: max iterations
% max_counter: tolerance for no improvements

if nargin == 2 || ~exist('max_counter','var')
    max_counter = 1000;
end

agr = agreement(c);
prob = agr/size(c,2);

ent = -nansum(prob.*log2(prob),2);
ent = ent.^-0.25;
entc = [0; cumsum(ent)];

klist = unique(max(c));
global_max = max(klist);

if length(unique(agr)) > 2
    k = mode(max(c));

    % initialize with a random partition from ensemble
    e_orig = c(:,randi(size(c,2)));
    [~,~,e] = unique(e_orig);  % relabel to contiguous integers
    [ari_vals,state] = init_ari_state(e,c,global_max);
    fitness = mean(ari_vals);

    n = length(e);
    internal_counter = 0;
    iter = 0;
    fitness_tracker = [0,fitness];
    while iter < numiter && internal_counter < max_counter
        iter = iter + 1;
        r = randi(n);              % pick a random node
        
        % r = sum(rand*entc(end) >= entc);


        a = agr(r,:);
        a(e == e(r)) = 0;          % can't map to own cluster
        if all(a==0), continue; end

        % sample a target node proportional to agreement
        p = [0,cumsum(a)];
        s = sum(rand*p(end) >= p);
        enew = e;

        oldC = e(r);
        newC = e(s);

        
        if oldC == newC, continue; end
        enew(r) = newC;

        % delta ARI update
        [ari_vals_new,state_new] = update_ari_state(r,oldC,newC,enew,c,ari_vals,state);
        fitnessnew = mean(ari_vals_new);

        internal_counter = internal_counter + 1;
        if any(length(unique(enew)) == klist)
            if (fitnessnew > fitness)
                % accept move
                fitness = fitnessnew;
                e = enew;
                ari_vals = ari_vals_new;
                state = state_new;
                internal_counter = 0;
                fitness_tracker = [fitness_tracker; iter,fitness];
            end
        end

        % % optional: monitor progress
        % if rand < 1/(numiter/1000)
        %     plot(iter,fitness,'ko'); hold on; drawnow;
        % end
    end
end
end

%% ---------- Helper functions ----------

function [ari_vals, state] = init_ari_state(e, C, global_max)
% initialize ARI state between candidate partition e and ensemble C
n = length(e);
m = size(C,2);

state = cell(m,1);
ari_vals = zeros(m,1);

for j = 1:m
    y = C(:,j);
    [ari_vals(j), state{j}] = init_one(e,y,global_max);
end
end

function [ari, st] = init_one(e,y,global_max)
n = length(e);

% Relabel to contiguous integers for accumarray
[~,~,eidx] = unique(e);
[~,~,yidx] = unique(y);

k = max(eidx);
l = max(yidx);

nij = accumarray([eidx,yidx],1,[global_max,global_max]);  % contingency table
ai = sum(nij,2);
bj = sum(nij,1);

st.nij = nij;
st.ai = ai;
st.bj = bj;
st.n = n;

ari = computeARI(nij,ai,bj,n);
end


function [ari_vals, state] = update_ari_state(r, oldC, newC, e, C, ari_vals, state)
m = size(C,2);

for j = 1:m
    y = C(r,j);   % cluster of r in ensemble partition j
    st = state{j};
    nij = st.nij; ai = st.ai; bj = st.bj; n = st.n;

    % decrement old
    nij(oldC,y) = nij(oldC,y) - 1;
    ai(oldC) = ai(oldC) - 1;

    % increment new
    nij(newC,y) = nij(newC,y) + 1;
    ai(newC) = ai(newC) + 1;

    ari_vals(j) = computeARI(nij,ai,bj,n);

    st.nij = nij; st.ai = ai; st.bj = bj;
    state{j} = st;
end
end

function ari = computeARI(nij,ai,bj,n)
nij2 = sum(sum(nij.^2 - nij))/2;
ai2 = sum(ai.^2 - ai)/2;
bj2 = sum(bj.^2 - bj)/2;
N2 = n*(n-1)/2;
expected = (ai2 * bj2) / N2;
maxden  = 0.5*(ai2 + bj2);
ari = (nij2 - expected) / (maxden - expected + eps);
end
