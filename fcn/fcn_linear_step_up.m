function [padj] = fcn_linear_step_up(p,q)
%FCN_LINEAR_STEP_UP         computes adjusted p-values
%
%   PADJ = FCN_LINEAR_STEP_SUP(P,Q) takes a vector of values, P, and an
%          accepted false discovery rate, Q, and computes an adjusted
%          p-value, PADJ.
%
%   INPUTS:     P,      vector of p values
%               Q,      false discovery rate
%
%   OUTPUTS: PADJ,      adjusted p value
%
%   Richard Betzel, Indiana University, 2012
%

%modification history
%02.18.2012 - original

p(isnan(p)) = [];
[r,c] = size(p);
if r > c
    p = p';
end

n = length(p);
t = (1:n)*(q/n);
p = sort(p,'ascend');

ind = find((t - p) > 0,1,'last');
padj = p(ind);

if isempty(padj)
    padj = 0;
end