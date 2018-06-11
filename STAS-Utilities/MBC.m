function [TpsiB,TBpsi] = MBC (N,b1,b2,b3,psi)
%
% Define the MBC transforms for a given trio of blade DOFs.
% Dimensions not listed in b1,b2,b3 are assigned the identity.
%
% Version:        Changes:
% --------        -------------
% 16.03.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 16.03.2018      Simple code, double-checked.
%
% Inputs:
% -------
% N               : Total dimension of the matrix.
% b1,b2,b3        : Vectors of length Np.
% psi             : Azimuth angle.
%
% Outputs:
% --------
% 

thrd = 1/3;
p2p3 = psi + 2*pi*thrd;
p4p3 = psi + 4*pi*thrd;
cp = [cos(psi);cos(p2p3);cos(p4p3)];
sp = [sin(psi);sin(p2p3);sin(p4p3)];
Bpsi = [thrd thrd thrd;2*thrd*cp.';2*thrd*sp.'];
psiB = [ones(3,1) cp sp];

Np = size(b1,1);

% And the "not blade" dimensions.
nb = ones(N,1);
nb(b1) = 0;
nb(b2) = 0;
nb(b3) = 0;
nb = logical(nb);
ind = [1:N].';
nind = ind(nb);
Nni = size(nind,1);

ii = [b1;b1;b1;b2;b2;b2;b3;b3;b3;nind];
jj = [b1;b2;b3;b1;b2;b3;b1;b2;b3;nind];
ss = [Bpsi(1,1)*ones(Np,1);Bpsi(1,2)*ones(Np,1);Bpsi(1,3)*ones(Np,1); ...
      Bpsi(2,1)*ones(Np,1);Bpsi(2,2)*ones(Np,1);Bpsi(2,3)*ones(Np,1); ...
      Bpsi(3,1)*ones(Np,1);Bpsi(3,2)*ones(Np,1);Bpsi(3,3)*ones(Np,1); ...
      ones(Nni,1)];
TBpsi = sparse(ii,jj,ss,N,N);
ss = [psiB(1,1)*ones(Np,1);psiB(1,2)*ones(Np,1);psiB(1,3)*ones(Np,1); ...
      psiB(2,1)*ones(Np,1);psiB(2,2)*ones(Np,1);psiB(2,3)*ones(Np,1); ...
      psiB(3,1)*ones(Np,1);psiB(3,2)*ones(Np,1);psiB(3,3)*ones(Np,1); ...
      ones(Nni,1)];
TpsiB = sparse(ii,jj,ss,N,N);
