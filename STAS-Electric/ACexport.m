function [dxdt,A,B] = ACexport (linFlag,x,u,p)
%
% State-space model of a generic AC transmission line with transformers
% and reactive compensation at each end.
%
% Bus 1 | Transformer | Bus 2 | Cable | Bus 3 | Transformer | Bus 4
%                    react.comp.     react.comp.
%
%   States:              y vector:             u vector:
%   i12d,q   1,2                               we          1
%   ir2d,q   3,4                               v1d,q      2,3
%   v2d,q    5,6                               v4d,q      4,5
%   i23d,q   7,8
%   ir3d,q   9,10
%   v3d,q   11,12
%   i34d,q  13,14
%   
%
% Version:        Changes:
% --------        -------------
% 08.03.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 08.03.2019      Derivatives verified by complex step.
%
% Inputs:
% -------
% p               :  1:  a12   (-)     Np/Ns turns ratio.
%                  2-5:  I12   (A)     RMS phase currents for nonlinear inductance
%                  6-9:  L12   (H)     Lp + (a^2)Ls nonlinear phase inductances
%                   10:  R12   (Ohms)  Phase resistance
%                   11:  L2    (H)     Effective inductance,  Bus 2 -> ground
%                   12:  R2    (Ohms)  Effective resistance,  Bus 2 -> ground
%                   13:  C2    (F)     Effective capacitance, Bus 2 -> ground
%                   14:  L23   (H)     Phase inductance in cable
%                   15:  R23   (Ohms)  Phase resistance in cable
%                   16:  L3    (H)     Effective inductance,  Bus 3 -> ground
%                   17:  R3    (Ohms)  Effective resistance,  Bus 3 -> ground
%                   18:  C3    (F)     Effective capacitance, Bus 3 -> ground
%                   19:  a34   (-)     Np/Ns turns ratio.
%                20-23:  I34   (A)     RMS phase currents for nonlinear inductance
%                24-27:  L34   (H)     Lp + (a^2)Ls nonlinear phase inductances
%                   28:  R34   (Ohms)  Phase resistance
%
% Outputs:
% --------
%

Nx = 14;
Nu = 5;

i12 = x(1:2);
ir2 = x(3:4);
v2  = x(5:6);
i23 = x(7:8);
ir3 = x(9:10);
v3  = x(11:12);
i34 = x(13:14);
we  = u(1);
v1  = u(2:3);
v4  = u(4:5);
a12 = p(1);
L2  = p(11);
R2  = p(12);
C2  = p(13);
L3  = p(16);
R3  = p(17);

x12 = i12;
u12 = [v1;v2;we];
p12 = p(1:10);
[dx12dt,y12out,aa12,bby12,cc12] = Transformer (linFlag,x12,u12,p12);

x23 = [i23;v3];
u23 = [i34+ir3;v2;we];
p23 = [p(15);p(14);p(18)];
[dx23dt,aa23,bby23] = PiLine (linFlag,x23,u23,p23);

x34 = i34;
u34 = [v3;v4;we];
p34 = p(19:28);
[dx34dt,y34out,aa34,bby34,cc34] = Transformer (linFlag,x34,u34,p34);

dxdt = zeros(Nx,1);
dxdt(1:2) = dx12dt;
dxdt([7 8 11 12]) = dx23dt;
dxdt(13:14) = dx34dt;

% Capacitance at Bus 2.
% C2 dv2/dt = -w C2 spin v2 + a12*i12 - ir2 - i23.
dxdt(5:6) = -we*[0 -1;1 0]*v2 + (a12*i12 - ir2 - i23)/C2;

% Inductance at Bus 2.
% L2 dir2/dt = - (w L2 spin + R) ir2 + v2.
dxdt(3:4) = -(we*[0 -1;1 0] + (R2/L2)*speye(2))*ir2 + v2/L2;

% Inductance at Bus 3.
% L3 dir3/dt = - (w L3 spin + R) ir3 + v3.
dxdt(9:10) = -(we*[0 -1;1 0] + (R3/L3)*speye(2))*ir3 + v3/L3;

if (linFlag == 1)

   A = spalloc(Nx,Nx,0.2*Nx*Nx);
   B = spalloc(Nx,Nu,0.1*Nx*Nu);

   ir = [1:2];
   A(ir,ir) = aa12;
   ic = [5:6];
   A(ir,ic) = bby12(:,3:4);
   ic = [2 3 1];
   B(ir,ic) = bby12(:,[1 2 5]);

   ir = [7 8 11 12];
   A(ir,ir) = aa23;
   ic = [13 14 9 10 5 6];
   A(ir,ic) = [bby23(:,1:2) bby23(:,1:2) bby23(:,3:4)];
   ic = 1;
   B(ir,ic) = bby23(:,5);

   ir = [13:14];
   A(ir,ir) = aa34;
   ic = [11:12];
   A(ir,ic) = bby34(:,1:2);
   ic = [4 5 1];
   B(ir,ic) = bby34(:,3:5);

   ir = [5:6];
   A(ir,ir) = -we*[0 -1;1 0];
   A(ir,1:2) = (a12/C2)*speye(2);
   A(ir,3:4) = (-1/C2)*speye(2);
   A(ir,7:8) = (-1/C2)*speye(2);
   B(ir,1) = -[0 -1;1 0]*v2;

   ir = [3:4];
   A(ir,ir) = -we*[0 -1;1 0] - (R2/L2)*speye(2);
   A(ir,5:6) = (1/L2)*speye(2);
   B(ir,1) = -[0 -1;1 0]*ir2;

   ir = [9:10];
   A(ir,ir) = -we*[0 -1;1 0] - (R3/L3)*speye(2);
   A(ir,11:12) = (1/L3)*speye(2);
   B(ir,1) = -[0 -1;1 0]*ir3;

else

   A = sparse(Nx,Nx);
   B = sparse(Nx,Nu);

end 

