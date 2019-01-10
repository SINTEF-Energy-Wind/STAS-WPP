function [dxdt,yout,A,By,C] = gridFrequency (Linflag,x,yin,params)
%
% Measurement of the grid frequency with a phase-locked loop.
%
%   states:           y vector:
%   th_m       1      th_e       1     in   (grid)
%   vmsd,q    2,3     vsd,q     2,3    in   (grid)
%   Psie       4      wem        4     out
%
% Version:        Changes:
% --------        -------------
% 25.11.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 25.11.2018      Derivatives verified with complex step.
%
% Inputs:
% -------
% Linflag         : Set to 1 for linearized equations.
% x               :  1    th_m  (rad)    Measured electrical angle.
%                   2,3   vms   (V)      Measured voltage at the network side terminals.
%                    4    Psie  (rad)    Integral of voltage angle.
% yin             :  1    th_e  (rad)    Actual grid voltage angle.
%                   2,3   vs    (V)      Voltage at the network side terminals.
% params          :  1    av    (rad/s)  Filter on voltage measurement.
%                    2    KPe   (rad/s)  Proportional gain.
%                    3    KIe   (rad/s2) Integral gain.
%                    4    weh   (rad/s)  Reference frequency.
%
% Outputs:
% --------
% yout            :  1    wem   (rad/s)  Measured electrical frequency.

Nx   = size(x,1);
Nyi  = size(yin,1);
Nyo  = 1;
Ny   = Nyi + Nyo;
dxdt = zeros(Nx,1);
A    = zeros(Nx,Nx);
By   = zeros(Nx,Ny);
C    = zeros(Ny,Nx);
Dy   = zeros(Ny,Ny);

av   = params(1);
KPe  = params(2);
KIe  = params(3);
weh  = params(4);

ce   = cos(yin(1));
se   = sin(yin(1));
cm   = cos(x(1));
sm   = sin(x(1));
cc   = ce*cm;
ss   = se*sm;
cs   = ce*sm;
sc   = se*cm;
Td   = [cc+ss, cs-sc; ...
        sc-cs, cc+ss];
vp   = Td*yin(2:3);

err       =  atan2c(x(3),x(2));
dxdt(2:3) = -av*x(2:3) + av*vp;
dxdt(4)   =  KIe*err;
yout      =  KPe*err + x(4) + weh;
dxdt(1)   =  yout;

if (Linflag == 1)

   dTde = [cs-sc, -cc-ss; ...
           cc+ss,  cs-sc];
   dTdm = [sc-cs,  cc+ss; ...
          -cc-ss,  sc-cs];

   nrm2 = x(2)^2 + x(3)^2;
   daty =  x(2)/nrm2;
   datx = -x(3)/nrm2;
   A  = [                    0, KPe*datx, KPe*daty, 1; ...
         av*dTdm(1,:)*yin(2:3),      -av,        0, 0; ...
         av*dTdm(2,:)*yin(2:3),        0,      -av, 0; ...
                             0, KIe*datx, KIe*daty, 0];
   By = [0, 0, 0, 0;                           ...
         av*dTde(1,:)*yin(2:3), av*Td(1,:), 0; ...
         av*dTde(2,:)*yin(2:3), av*Td(2,:), 0; ...
         0, 0, 0, 0];
   C  = [zeros(3,4); ...
         0, KPe*datx, KPe*daty, 1];

end
