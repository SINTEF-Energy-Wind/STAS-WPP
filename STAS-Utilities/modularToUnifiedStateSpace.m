function [A,B,C,D] = modularToUnifiedStateSpace(Am,Bu,By,Cm,Du,Dy, ...
                                                Dnorm)
%
% Converts a modular state space of the form
% dx/dt = Am x + Bu u + By y
%   y   = Cm x + Du u + Dy y
% into the standard form
% dx/dt = Ax + Bu
%   y   = Cx + Du
%
% Dnorm is a vector of length Ny which is used to balance Dy, to
% avoid some numerical conditioning warnings.  If unused, set this
% to ones(Ny,1).
%

Nx = size(Am,1);
Ny = size(Cm,1);
Nu = size(Bu,2);

A = sparse(Nx,Nx);
B = sparse(Nx,Nu);
C = sparse(Ny,Nx);
D = sparse(Ny,Nu);

ii = [1:Ny];
jj = [1:Ny];
ImDy = sparse(ii,jj,ones(Ny,1)) - Dy;

% inv(ImDy) may be such that there are warnings of ill-conditioning, 
% however in all cases so far the inverse is precise, as verified by
% ImDy1*ImDy and ImDy*ImDy1.  That said, here is a way to avoid the
% warnings.
DD = sparse(ii,jj,Dnorm);
iDD = inv(DD);
ImDy1 = DD*((iDD*ImDy*DD)\iDD);

C(:,:) = ImDy1*Cm;
D(:,:) = ImDy1*Du;

A(:,:) = Am + By*C;
B(:,:) = Bu + By*D;

%del = 1e-12;
%A(abs(A)<del) = 0;
%B(abs(B)<del) = 0;
%C(abs(C)<del) = 0;
%D(abs(D)<del) = 0;
