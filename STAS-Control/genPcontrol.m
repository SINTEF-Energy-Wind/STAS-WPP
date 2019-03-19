function [dxdt,yout,A,B,C,D] = genPcontrol (x,u,params)
%
% Feedback control on active power, taking the measured and commanded
% power as inputs, and outputting a q-axis generator current command
% to be fed to the electrical controls.
%
%   states:             y vector:           u vector:
%   Pem     1           igq     1           Phat    1  (Turbine control)
%   PsiP    2                               Pe      2  (Electrical power)
%   nd     3,4
%    
% Version:        Changes:
% --------        -------------
% 25.01.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 25.01.2019      Derivatives verified by complex step.
%
% Inputs:
% -------
% params          : 1: Kp    (A/W)      Proportional gain.
%                   2: Ki    (A/Ws)     Integral gain.
%                   3: ap    (rad/s)    LP filter on power.
%                   4: anp   (rad/s)    Notch filter frequency.
%                   5: z1np  (-)        Notch filter numerator.
%                   6: z2np  (-)        Notch filter denominator.
%
% Outputs:
% --------
% State matrices.

Pem  = x(1);
PsiP = x(2);
nd   = x(3:4);
Phat = u(1);
Pe   = u(2);
Kp   = params(1);
Ki   = params(2);
ap   = params(3);
anp  = params(4);
z1np = params(5);
z2np = params(6);

tz2a = 2*z2np*anp;
tzza = 2*(z1np - z2np)*anp;
a2   = anp^2;

dxdt = zeros(4,1);
dxdt(1) = -ap*Pem + ap*Pe;

epsP = Phat - Pem;
dxdt(3) = nd(2);
dxdt(4) = -a2*nd(1) - tz2a*nd(2) + epsP;

yfilt = tzza*nd(2) + epsP;

dxdt(2) = Ki*yfilt;

yout = Kp*yfilt + PsiP;

A  = [-ap, 0, 0, 0;    ...
      -Ki, 0, 0, Ki*tzza; ...
       0, 0, 0, 1;     ...
      -1, 0, -a2, -tz2a];
B  = [0, ap; ...
      Ki, 0; ...
      0, 0;  ...
      1, 0];
C  = [-Kp, 1, 0, Kp*tzza];
D  = [Kp, 0];

