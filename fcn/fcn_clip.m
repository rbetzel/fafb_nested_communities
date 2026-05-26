function x = fcn_clip(x,lims)
x(x < lims(1)) = lims(1);
x(x > lims(2)) = lims(2);