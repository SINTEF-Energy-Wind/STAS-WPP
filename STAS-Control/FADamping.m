function [dxdt,yout,A,B,C,D] = FADamping (x,u,p,KFTab)
%
% This function builds a fore-aft tower damping control loop.  This
% is nominally based on an accelerometer measurement in the nacelle,
% but for computational reasons (avoiding inverting the mass matrix)
% this implementation takes the velocity directly from the
% structural model.
%
%   States:              y vector:             u vector:
%   vmF    1             bF       1            vF      1
%   ivmF   2                                   betm    2
%   vdF    3
%
% The input betm is used only as a scheduling variable.  It is named 
% as if it is the blade pitch, but it doesn't need to be.  It just
% needs to be defined consistently with the scheduling table KFTab.
% An alternative scheduling variable would be the windspeed.
%
% Version:        Changes:
% --------        -------------
% 02.02.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 02.02.2019      Derivatives verified with complex step.  (Using the
%                 form of gains1 that is suitable for complex step.)
%
% Inputs:
% -------
% x               : 1: vmF   (m/s)     Measured nacelle FA velocity.
%                   2: ivmF  (m)       Integral of vmF.
%                   3: vdF   (m/s)     Bandpass filtered velocity.
% u               : 1: vF    (m/s)     Nacelle FA velocity.
%                   2: betm  (rad)     blade pitch angle.
% p               : 1: aF    (rad/s)   Frequency in rad/s.
%                   2: zetaF (-)       Damping ratio.
%                   3: azF   (rad/s)   Numerator frequency of phase shift.
%                   4: apF   (rad/s)   Denominator frequency of phase shift.
% KFTab           : Table of KF(beta), scheduled gains.

Nx = 3;
Nu = 2;
Ny = 1;

vmF   = x(1);
ivmF  = x(2);
vdF   = x(3);
vF    = u(1);
betm  = u(2);
aF    = p(1);
zetaF = p(2);
azF   = p(3);
apF   = p(4);

[KF,dKF] = gains1 (betm,KFTab);

A  = sparse(Nx,Nx);
B  = sparse(Nx,Nu);
C  = sparse(Ny,Nx);
D  = sparse(Ny,Nu);

tza = 2*zetaF*aF;
apaz = apF/azF;
a2 = aF^2;

dxdt = [-2*zetaF*aF*vmF - (aF^2)*ivmF + 2*zetaF*aF*vF; ...
        vmF;                                           ...
        (-tza*apaz + apF)*vmF - apaz*a2*ivmF - apF*vdF + tza*apaz*vF];
yout = KF*vdF;

% divmF/dt = vmF.
A(2,1) = 1;

% dvmF/dt = -2*zetaF*aF*vmF - (aF^2)*ivmF + 2*zetaF*aF*vF.
A(1,1) = -tza;
A(1,2) = -a2;
B(1,1) =  tza;

% dvdF/dt = (-tza*apaz + ap)*vmF - apaz*a2*ivmF - ap*vdF
%         + tza*apaz*vF.
A(3,1) = -tza*apaz + apF;
A(3,2) = -apaz*a2;
A(3,3) = -apF;
B(3,1) =  tza*apaz;

% beta = KF*vdF, since +beta means -Ft, and we want Ft to go down
% if vF is positive (moving downwind).
C(1,3) = KF;
D(1,2) = dKF*vdF;

% db/dt = KF*dvdF/dt.
%C(2,:) = A(3,:)*KF;
%D(2,1) = B(3,1)*KF;
%D(2,2) = dKF*dxdt(3);


