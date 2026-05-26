function lims = fcn_axlims(x,y,offset)
if nargin == 2
    offset = 0.05;
end
x = x(:); y = y(:);
x(isnan(x) | isinf(x)) = [];
y(isnan(y) | isinf(y)) = [];
set(gca,...
    'xlim',[[min(x),max(x)] + offset*[-1,1]*range(x)],...
    'ylim',[[min(y),max(y)] + offset*[-1,1]*range(y)]);

lims(1,:) = get(gca,'xlim');
lims(2,:) = get(gca,'ylim');