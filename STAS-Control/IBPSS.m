function [dxdt,yout,A,B,C,D] = IBPSS (x,u,p,KSTab)
%
% This function builds an IBP pitch controller which acts to damp
% side-to-side motion of the tower.  The sine (q) component of pitch
% acts to provide a rotor in-plane force component that opposes the
% tower motion over a relevant frequency band.
%
%   States:              y vector:             u vector:
%   vmS      1           bet      1:3          vS     1
%   ivmS     2                                 azi    2
%   vdS      3                                 sch    3
%                                              
%                       
% Version:        Changes:
% --------        -------------
% 04.02.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 04.02.2019      Derivatives verified by complex step.
%
% Inputs:
% -------
% x               : 1: vmS    (m/s)    Velocity of tower side-to-side motion.
%                   2: ivmS   (m)      Integral of vmS.
%                   3: vdS    (m/s)    Bandpass filtered velocity.
% u               : 1: vS    (m/s)     Nacelle SS velocity.
%                   2: sch             Gain scheduling variable.
% p               : 1: aS    (rad/s)   Frequency in rad/s.
%                   2: zetaS (-)       Damping ratio.
%                   3: azS   (rad/s)   Numerator frequency of phase shift.
%                   4: apS   (rad/s)   Denominator frequency of phase shift.
% KSTab           : Table of KS(sch), scheduled gains.

Nx = 3;
Nu = 3;
Ny = 3;

vmS   = x(1);
ivmS  = x(2);
vdS   = x(3);
vS    = u(1);
azi   = u(2);
sch   = u(3);
aS    = p(1);
zetaS = p(2);
azS   = p(3);
apS   = p(4);

[KS,dKS] = gains1 (sch,KSTab);

A  = sparse(Nx,Nx);
B  = sparse(Nx,Nu);
C  = sparse(Ny,Nx);
D  = sparse(Ny,Nu);

[TpsiB,TBpsi] = MBC (3,1,2,3,azi);
[dTpsiB,dTBpsi] = derivMBC (3,1,2,3,azi);

tza = 2*zetaS*aS;
apaz = apS/azS;
a2 = aS^2;

dxdt = [-2*zetaS*aS*vmS - (aS^2)*ivmS + 2*zetaS*aS*vS; ...
        vmS;                                           ...
        (-tza*apaz + apS)*vmS - apaz*a2*ivmS - apS*vdS + tza*apaz*vS];
yout = -KS*vdS*TpsiB(:,3);

% divmS/dt = vmS.
A(2,1) = 1;

% dvmS/dt = -2*zetaS*aS*vmS - (aS^2)*ivmS + 2*zetaS*aS*vS.
A(1,1) = -tza;
A(1,2) = -a2;
B(1,1) =  tza;

% dvdS/dt = (-tza*apaz + ap)*vmS - apaz*a2*ivmS - ap*vdS
%         + tza*apaz*vS.
A(3,1) = -tza*apaz + apS;
A(3,2) = -apaz*a2;
A(3,3) = -apS;
B(3,1) =  tza*apaz;

% TgS = -KS*vdS*(T_psi^B)q.
C(:,3) = -KS*TpsiB(:,3);
D(:,2) = -KS*vdS*dTpsiB(:,3);
D(:,3) = -dKS*vdS*TpsiB(:,3);



