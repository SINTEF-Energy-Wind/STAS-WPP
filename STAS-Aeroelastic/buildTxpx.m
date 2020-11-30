function [x,T,dT] = buildTxpx (xpsi,dxpsi,b1,b2,b3,iazi)
%
% Transform of state variables from MBC to body coordinates.  Also
% computes the values of the states in body coordinates.
%
% By definition, [eta;etad] = T_psi^B [eta^psi;etad^psi].  Then,
% dx/dt = T_xp^x dx^psi/dt.
%
% Version:        Changes:
% --------        -------------
% 14.01.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 14.01.2020      Derivatives of Txpx verified by complex step.
%
% Inputs:
% -------
% xpsi            : State vector in MBC coordinates.
% dxpsi           : d xpsi/dt, only used for dT, not for T.
% b1,b2,b3        : DOFs corresponding to blades 1,2,3.
% psi0            : Reference azimuth angle for transform (should be
%                   independent of this).
% iazi            : The index of the rotor azimuth DOF in x^psi.
%
% Outputs:
% --------
% x               : x in body coordinates.
% T               : T_xpsi^x.
% dT              : = (dT_xpsi^x)*dxpsi. 

azi = xpsi(iazi);

Nx = size(xpsi,1);
[TpsiB,TBpsi]     = MBC (Nx,b1,b2,b3,azi);
[dTpsiB,dTBpsi]   = derivMBC (Nx,b1,b2,b3,azi);
[d2TpsiB,d2TBpsi] = secondDerivMBC (Nx,b1,b2,b3,azi);

x = TpsiB*xpsi;

dTx = dTpsiB*xpsi;
T = TpsiB;
T(:,iazi) = T(:,iazi) + dTx;

W   = dxpsi(iazi);
dT  = W*dTpsiB;
dT(:,iazi) = dT(:,iazi) + dTpsiB*dxpsi + W*d2TpsiB*xpsi;

