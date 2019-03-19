function [dxdt,yout,A,B,C,D] = SSDamping (x,u,p,KSTab)
%
% This function builds a side-to-side tower damping control loop.
% The acceleration in the side-to-side direction is measured at the
% nacelle.  This is integrated to get the velocity, then band-pass
% filtered at the tower frequency.
%
% The present implementation takes nacelle velocity as input, for
% computational reasons (avoiding inverting the mass matrix).
%
%   States:              y vector:             u vector:
%   vmS      1           TgS      1            vS       1
%   ivmS     2                                 sch      2
%   vdS      3
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
% x               : 1: vmS   (m/s)     Measured nacelle SS velocity.
%                   2: ivmS  (m)       Integral of vmS.
%                   3: vdS   (m/s)     Bandpass filtered velocity.
% u               : 1: vS    (m/s)     Nacelle SS velocity.
%                   2: sch             Gain scheduling variable.
% p               : 1: aS    (rad/s)   Frequency in rad/s.
%                   2: zetaS (-)       Damping ratio.
%                   3: azS   (rad/s)   Numerator frequency of phase shift.
%                   4: apS   (rad/s)   Denominator frequency of phase shift.
% KSTab           : Table of KS(sch), scheduled gains.

Nx = 3;
Nu = 2;
Ny = 1;

vmS   = x(1);
ivmS  = x(2);
vdS   = x(3);
vS    = u(1);
sch   = u(2);
aS    = p(1);
zetaS = p(2);
azS   = p(3);
apS   = p(4);

[KS,dKS] = gains1 (sch,KSTab);

A  = sparse(Nx,Nx);
B  = sparse(Nx,Nu);
C  = sparse(Ny,Nx);
D  = sparse(Ny,Nu);

tza = 2*zetaS*aS;
apaz = apS/azS;
a2 = aS^2;

dxdt = [-2*zetaS*aS*vmS - (aS^2)*ivmS + 2*zetaS*aS*vS; ...
        vmS;                                           ...
        (-tza*apaz + apS)*vmS - apaz*a2*ivmS - apS*vdS + tza*apaz*vS];
yout = KS*vdS;

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

% TgS = KS*vdS.
C(1,3) = KS;
D(1,2) = dKS*vdS;

