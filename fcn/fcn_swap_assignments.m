function [A,M] = fcn_swap_assignments(C,A,M,i,old,new)
Ci = C(:,i);
A(i,old) = 0;
A(i,new) = 1;
M(:,old) = M(:,old) - Ci;
M(:,new) = M(:,new) + Ci;