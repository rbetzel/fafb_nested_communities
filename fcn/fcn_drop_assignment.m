function [A,M] = fcn_drop_assignment(C,A,M,i,old)
Ci = C(:,i);
A(i,old) = 0;
M(:,old) = M(:,old) - Ci;