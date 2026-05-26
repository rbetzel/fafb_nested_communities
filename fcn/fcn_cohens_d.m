function d = fcn_cohens_d(data, group)
% COHENS_D Compute Cohen's d for two groups
%
%   d = COHENS_D(data, group)
%
%   INPUTS:
%       data  - numeric vector of observations
%       group - vector of group labels (must contain exactly two unique values)
%
%   OUTPUT:
%       d     - Cohen's d (standardized mean difference)
%
%   EXAMPLE:
%       data = [2.1 2.5 2.9 3.0 3.1 5.2 5.5 5.8 6.0 6.1];
%       group = [1 1 1 1 1 2 2 2 2 2];
%       d = cohens_d(data, group)

% Check input lengths
if numel(data) ~= numel(group)
    error('data and group must have the same length.');
end

% Identify groups
groups = unique(group);
if numel(groups) ~= 2
    error('group must contain exactly two unique values.');
end

% Split data into two groups
x1 = data(group == groups(1));
x2 = data(group == groups(2));

% Sample sizes
n1 = numel(x1);
n2 = numel(x2);

% Means and standard deviations
m1 = mean(x1);
m2 = mean(x2);
s1 = std(x1, 1); % population std (optional)
s2 = std(x2, 1);

% Use sample SD (set flag 0 instead of 1)
s1 = std(x1, 0);
s2 = std(x2, 0);

% Pooled standard deviation
sp = sqrt(((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2));

% Cohen's d
d = (m1 - m2) / sp;
