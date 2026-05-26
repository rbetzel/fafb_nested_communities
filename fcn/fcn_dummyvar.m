function I = fcn_dummyvar(val,num_cols,normflag)
% dummy code categorical data

unq = unique(nonzeros(val));
if ~exist('num_cols','var') | isempty(num_cols)
    I = zeros(size(val,1),max(unq));
else
    I = zeros(size(val,1),num_cols);
end
for i = 1:length(unq)
    I(:,unq(i)) = val == unq(i);
end
if ~exist('normflag','var')
    normflag = false;
end
if normflag
    I = bsxfun(@rdivide,I,sum(I));
end
