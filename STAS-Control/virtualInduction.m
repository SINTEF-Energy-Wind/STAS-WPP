function [dxdt,yout,A,B,C] = virtualInduction (x,u,p)
%
% This function implements a virtual induction generator for
% damping driveshaft torsional resonance.  The electrical speed
% is band-pass filtered to produce an in-phase torque signal,
% around the resonant frequency band.
%
%   States:              y vector:             u vector:
%   wgbar   1            Tgi    1              wg     1
%   iwgb    2
%
% Version:        Changes:
% --------        -------------
% 02.02.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 02.02.2019      Verified derivatives with complex step.  Verified
%                 transfer function.
%
% Inputs:
% -------
% x               : 1: wgbar  (rad/s)   filtered gen. speed.
%                   2: iwgb   (rad)     integrated wgbar.
% u               : 1: wg     (rad/s)   gen. rotor speed.
% p               : 1: ag     (rad/s)   BP filter frequency.
%                   2: zetag  (-)       BP filter damping.
%                   3: Kd     (Nms/rad) Induction generator stiffness.

Nx = 2;
Nu = 1;
Ny = 1;

ag    = p(1);
zetag = p(2);
Kd    = p(3);

A  = sparse(Nx,Nx);
B  = sparse(Nx,Nu);
C  = sparse(Ny,Nx);

% diwgbar/dt = wgbar.
A(2,1) = 1;

% dwgbar/dt = -2*zetag*ag*wgbar - (ag^2)*iwgb + 2*zetag*ag*wg
A(1,1) = -2*zetag*ag;
A(1,2) = -(ag^2);
B(1,1) =  2*zetag*ag;

% Tgi = Kd*wgbar.
C(1,1) = Kd;

% The equations are linear.
dxdt = A*x + B*u;
yout = C*x;
