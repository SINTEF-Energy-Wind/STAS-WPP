function Am1B = balanceDiv (A,B)
%
% Compute a balanced solution to A\B, handy when A is not well
% conditioned.
%
% To compute an inverse, enter B = eye(size(A)).
%

Nx = size(A,1);

[Diag,Perm,Anorm] = balance(A,'noperm');
ii = [1:Nx];
jj = ii;
DD = sparse(ii,jj,Diag,Nx,Nx);
iDD = inv(DD);

%cond(Anorm)

Am1B = DD*(Anorm\(iDD*B));
