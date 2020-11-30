function [A,B] = Vcascade (NvC)

ii = [2:NvC];
jj = [1:NvC-1];
ss = ones(NvC-1,1);
A = sparse(ii,jj,ss,NvC,NvC);
B = sparse(NvC,1);
B(1) = 1;
