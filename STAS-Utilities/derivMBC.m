function [dTpsiB,dTBpsi] = derivMBC (N,b1,b2,b3,psi)
%
% Define the derivative of the MBC transforms for a given trio
% of blade DOFs.
%
% Version:        Changes:
% --------        -------------
% 20.04.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 20.04.2018      Identical to MBC, except for derivatives dBpsi
%                 and dpsiB.
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
dBpsi = [zeros(1,3);-2*thrd*sp.';2*thrd*cp.'];
dpsiB = [zeros(3,1) -sp cp];

Np = size(b1,1);

ii = [b1;b1;b1;b2;b2;b2;b3;b3;b3];
jj = [b1;b2;b3;b1;b2;b3;b1;b2;b3];
ss = [dBpsi(1,1)*ones(Np,1);dBpsi(1,2)*ones(Np,1);dBpsi(1,3)*ones(Np,1); ...
      dBpsi(2,1)*ones(Np,1);dBpsi(2,2)*ones(Np,1);dBpsi(2,3)*ones(Np,1); ...
      dBpsi(3,1)*ones(Np,1);dBpsi(3,2)*ones(Np,1);dBpsi(3,3)*ones(Np,1)];
dTBpsi = sparse(ii,jj,ss,N,N);
ss = [dpsiB(1,1)*ones(Np,1);dpsiB(1,2)*ones(Np,1);dpsiB(1,3)*ones(Np,1); ...
      dpsiB(2,1)*ones(Np,1);dpsiB(2,2)*ones(Np,1);dpsiB(2,3)*ones(Np,1); ...
      dpsiB(3,1)*ones(Np,1);dpsiB(3,2)*ones(Np,1);dpsiB(3,3)*ones(Np,1)];
dTpsiB = sparse(ii,jj,ss,N,N);
