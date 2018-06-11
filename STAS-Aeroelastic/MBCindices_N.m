function [b1,b2,b3] = MBCindices_N (N)

b1 = [1:N].';
b2 = N + [1:N].';
b3 = 2*N + [1:N].';
