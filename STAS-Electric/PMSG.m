function [dxdt,yout,LA,LBy,C] = PMSG (Linflag,x,yin,params)
%
% The generator is modelled with an equivalent electric circuit
% in the dq frame.  It is assumed that the inductance and
% resistance matrices in the abc frame are diagonal, that is,
% no magnetic coupling between the phases.
%
%   states:           y vector:
%   igd,q     1,2     vgd,q     1,2    in   (gen. current control)
%                     wg         3     in   (shaft)
%                     Tg         4     out
%
% ig is the stator current, vg is the terminal voltage, wg is the
% electrical speed, and Tg is the mechanical torque.
%
% Version:        Changes:
% --------        -------------
% 02.05.2016      Code based on buildGenerator.m.
% 11.07.2017      Adapted for complex step derivatives.
% 27.08.2018      Modified for linear/nonlinear equation pairs
%                 and nonlinear inductance.
%
% Version:        Verification:
% --------        -------------
% 02.05.2016      
% 11.07.2017      
% 27.08.2018      Derivatives verified by complex step.
%
% Inputs:
% -------
% Linflag         : Set to 1 for linearized equations.
% x               : 1,2:  ig    (A)     Stator currents
% yin             : 1,2:  vg    (V)     Generator terminal voltage
%                   3:    wg    (rad/s) Generator electrical speed
% params          : 1:    np    (-)     Number of poles.
%                   2-5:  Ig    (A)     RMS phase currents for nonlinear inductance.
%                   6-9:  Lg    (H)     Nonlinear phase inductances.
%                   10:   Rg    (Ohms)  Phase resistance.
%                   11:   lamr  (Wb)    d-component flux linkage per phase.
%
% Outputs:
% --------
% yout            : 1    Tg    (Nm)    Generator torque.

Nx = 2;
Ny = 4;

dxdt = zeros(Nx,1);
yout = 0;

L = sparse(Nx,Nx);
A = sparse(Nx,Nx);
By = sparse(Nx,Ny);
C = sparse(Ny,Nx);

ig = x;

vg = yin(1:2);
wg = yin(3);

np  = params(1);
IIg = params(2:5);
Lg  = params(6:9);
Rg  = params(10);
lam = params(11);

two = 2*pi/3;
four = 4*pi/3;
sq = sqrt(2/3);
sq2 = sq/sqrt(2);

th = 0;

ct = cos(th);
st = sin(th);
c2 = cos(th-two);
s2 = sin(th-two);
c4 = cos(th-four);
s4 = sin(th-four);

lamr = lam*[ct c2 c4].';

Tat = sq*[ct  c2  c4; ...
         -st -s2 -s4];
Tta = sq*[ct -st; ...
          c2 -s2; ...
          c4 -s4];
dTtadth = sq*[-st -ct; ...
              -s2 -c2; ...
              -s4 -c4];

% Compute the RMS phase current and phase inductance.
magI = sqrt(ig(1)^2 + ig(2)^2);
angi = atan2c(ig(2),ig(1));
im1  = cos(angi);
im2  = sin(angi);
Irms = magI*sq2;
[Lph,dLph] = inductance (IIg,Lg,Irms);

% dq-transformed permanent magnet flux linkage.
lamdq   = Tat*lamr;
TTL     = Tat*dTtadth*lamdq;
WTL     = wg*TTL;
Lambda  = Lph*eye(3);
dLambda = dLph*eye(3);
dLamN   = [dLambda*sq2*im1 dLambda*sq2*im2];
TLT     = Tat*Lambda*Tta;
TLamDT  = Tat*Lambda*dTtadth;
WTLamDT = wg*TLamDT;
RR      = Rg*eye(3);
TRT     = Tat*RR*Tta;

EW = -0.5*np*TTL;

dxdt = -TLT\((TRT + WTLamDT)*ig + vg + WTL);
yout = (EW.')*ig;

if (Linflag == 1)

   L(:,:) = Tat*Lambda*Tta;

   A(:,:) = -(TRT + WTLamDT);
   vec = Tta*dxdt + wg*dTtadth*ig;
   A(:,1) = A(:,1) - Tat*dLamN(:,1:3)*vec;
   A(:,2) = A(:,2) - Tat*dLamN(:,4:6)*vec;

   By(:,1:2) = -eye(2);
   By(:,3) = -TLamDT*ig - TTL;
   C(4,:) = EW.';

   LA = L\A;
   LBy = L\By;

else

   LA = sparse(Nx,Nx);
   LBy = sparse(Nx,Ny);

end

%{
Tat
lam
lamr
wg
TLT
WTLamDT
TRT
vg
WTL
lamdq
%}