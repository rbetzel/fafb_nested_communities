function sz = fcn_sz(k,szrng)
sz = interp1([min(k),max(k)],szrng,k);