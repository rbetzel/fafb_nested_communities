function cmap = fcn_cmap(base,ndiv)
if nargin == 1
    ndiv = 256;
end
cmap = interp1(linspace(0,1,size(base,1)),base,linspace(0,1,ndiv));