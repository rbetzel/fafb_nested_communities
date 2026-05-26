function [xx,yy] = fcn_get_boxplot_skeleton(x,lab)

w = 1.5;

iqr = [25,50,75];
wid = 0.3;

u = unique(lab);

xx = cell(length(u),1);
yy = xx;
for j = 1:length(u)
    
    idx = lab == u(j);
    vals = x(idx);
    prct = prctile(vals,iqr);
    
    dr = prct(3) - prct(1);
    upper = prct(3) + w*dr;
    lower = prct(1) - w*dr;
    
    whisk1 = min(vals(vals >= lower));
    whisk2 = max(vals(vals <= upper));

    xvals1 = j + [-wid,+wid,+wid,-wid,-wid];
    yvals1 = [prct(1),prct(1),prct(3),prct(3),prct(1)];

    xvals2 = j + [-wid,wid];
    yvals2 = prct(2)*ones(1,2);

    xvals3 = [j*ones(1,2), nan, j*ones(1,2)];
    yvals3 = [prct(1),whisk1,nan,prct(3),whisk2];

    xvals = [xvals1,nan,xvals2,nan,xvals3];
    yvals = [yvals1,nan,yvals2,nan,yvals3];

    xx{j} = xvals;
    yy{j} = yvals;

end