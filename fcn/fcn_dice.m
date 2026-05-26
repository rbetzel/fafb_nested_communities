function d = fcn_dice(M,L)
d = 2*sum(M & L)./(sum(M) + sum(L));