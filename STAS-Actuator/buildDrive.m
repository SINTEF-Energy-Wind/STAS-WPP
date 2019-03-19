function [dxdt,yout,A,By,C,Dy] = buildDrive (Linflag,x,yin,params)
%
% This builds the equations associated with a pitch or yaw drive.
%
%   states:             y vector:
%   ba      1           bhat      1          in   (Controller)
%   dbadt   2           b,dbdt   2,3         in   (Turbine)
%                       Ta        4         out
%    
% Version:        Changes:
% --------        -------------
% 11.01.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 11.01.2019      Derivatives verified by complex step.
%
% Inputs:
% -------
% params          : 1: ka    (Nm/rad)   Torsional stiffness.
%                   2: ca    (Nms/rad)  Torsional damping.
%                   3: aa    (rad/s)    Frequency of motor dynamics.
%                   4: K     (1/s)      Gain on pitch rate control.
%                   5: bmax  (rad)      Max slew angle.
%                   6: bmin  (rad)      Min slew angle.
%                   7: bdmax (rad/s)    Max slew rate.
%
% Outputs:
% --------
% State matrices.

Nx = 2;
Ny = 4;

ka    = params(1);
ca    = params(2);
aa    = params(3);
K     = params(4);
bmax  = params(5);
bmin  = params(6);
bdmax = params(7);

h = 0.5*(bmax - bmin);
g = bmin + h;

[S1,dS1] = saturate ((yin(1)-g)/h,[0.9;1.1]);
bst = h*S1 + g;
epsd = bst - x(1);
badh = K*epsd;
[S2,dS2] = saturate (badh/bdmax,[0.8,1.2]);
badhs = bdmax*S2;
dxdt = [x(2);(-aa*x(2) + aa*badhs)];
yout = ka*(x(1) - yin(2)) + ca*(x(2) - yin(3));

adSK = aa*dS2*K;

if (Linflag == 1)

   A = [0, 1; -adSK, -aa];
   By = [[0; adSK*dS1], zeros(2,3)];
   C = [zeros(3,2);[ka, ca]];
   Dy = [zeros(3,4);[0, -ka, -ca, 0]];

else

   A = zeros(2,2);
   By = zeros(2,4);
   C = zeros(4,2);
   Dy = zeros(4,4);

end

