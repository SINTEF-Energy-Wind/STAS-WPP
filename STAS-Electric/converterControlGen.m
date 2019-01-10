function [dxdt,yout,A,By,C,Dy] = converterControlGen (Linflag,x,yin,params)
%
% The generator-side converter is controlled so as to provide a commanded
% active power, and zero d-axis current.  With the present convention
% used for the generator, with the d axis following the rotor magnetic
% N pole and the assumed direction of the windings, a positive (braking)
% torque is obtained with a negative iq.
%
%   states:           y vector:
%   imgd,q    1,2     wg         1     in   (shaft)
%   Psig      3:4     igd,q     2,3    in   (generator)
%   wemg       5      ihgd,q    4,5    in   (active power control)
%                     vgd,q     6,7    out
%
% Version:        Changes:
% --------        -------------
% 02.05.2016      Code based on buildGenerator.m.
% 11.07.2017      Adapted for complex step derivatives.
% 27.08.2018      Modified for linear/nonlinear equation pairs and
%                 nonlinear inductance.
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
% x               : 1,2:  img   (A)     Measured stator currents.
%                   3,4:  Psig  (A/s)   Integrated current error KI*eps.
%                   5:    wemg  (rad/s) Measured generator electrical speed.
% yin             : 1:    wg    (rad/s) Generator electrical speed.
%                   2,3:  ig    (A)     Generator stator current.
%                   4,5:  ihg   (A)     Current command from act. pow. control.
% params          : 1-4:  Ig    (A)     RMS phase currents for nonlinear inductance.
%                   5-8:  Lg    (H)     Nonlinear phase inductances.
%                   9:    lamr  (Wb)    d-component flux linkage per phase.
%                   10:   KP    (1/s)   Proportional gain.
%                   11:   KI    (1/s^2) Integral gain.
%                   12:   KF    (-)     Feed-forward gain on flux linkage.
%                   13:   ag    (rad/s) Generator current filter frequency.
%                   14:   aw    (rad/s) Generator electric speed filter frequency.
%
% Outputs:
% --------
% yout            : 1,2   vg    (V)     Generator terminal voltage.

Nx  = 5;
Ny  = 7;
Nyi = 5;
Nyo = 2;

dxdt = zeros(Nx,1);
yout = zeros(2,1);
C    = zeros(Ny,Nx);
Dy   = zeros(Ny,Ny);

I0   = params(1:4);
L0   = params(5:8);
lamr = params(9);
KP   = params(10);
KI   = params(11);
KF   = params(12);
ag   = params(13);
aw   = params(14);

A    = [-ag 0  0  0  0; ...
         0 -ag 0  0  0; ...
        -KI 0  0  0  0; ...
         0 -KI 0  0  0; ...
         0  0  0  0 -aw];
By   = [0  ag 0  0  0  0  0; ...
        0  0  ag 0  0  0  0; ...
        0  0  0  KI 0  0  0; ...
        0  0  0  0  KI 0  0; ...
        aw 0  0  0  0  0  0];

dxdt = A*x + By(:,1:Nyi)*yin;

epsg = yin(4:5) - x(1:2);
IIm = sqrt(x(1)^2 + x(2)^2);
[Lm,dLm] = inductance (I0,L0,IIm);
term = KP*epsg + x(3:4) + x(5)*[0,-1;1,0]*x(1:2) + KF*x(5)*lamr;
yout = -Lm*term;

if (Linflag == 1)

   C(Nyi+[1:2],1:2)  =  Lm*([KP 0;0 KP] - x(5)*[0,-1;1,0]) ...
                     -  term*dLm*[x(1) x(2)]/IIm;
   C(Nyi+[1:2],3:4)  = -Lm*eye(2);
   C(Nyi+[1:2],5)    = -Lm*([0,-1;1,0]*x(1:2) + KF*lamr);

   Dy(Nyi+[1:2],4:5) = -Lm*[KP 0;0 KP];

end

