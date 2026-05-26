function [cinew,idx] = fcn_sort_communities(ci)
% relabel communities based on size
%
%   Richard Betzel, 2025, University of Minnesota
%
h = hist(ci,1:max(ci));
[~,idx] = sort(h,'descend');
cinew = zeros(size(ci));
for j = 1:max(ci)
    jdx = ci == idx(j);
    cinew(jdx) = j;
end