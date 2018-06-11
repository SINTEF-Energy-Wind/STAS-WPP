function [Viq,C,Dy] = Vi_qs (Nb,Viz,Vzts,rho,r,L,f,Fzts)
%
% Quasi-steady induced velocity.
%
%   States:           y vector:         u vector:
%   Viz               Vzts   1:3
%                     r       4
%                     L       5
%                     f       6
%                     Fzts   7:9
%
% Version:        Changes:
% --------        -------------
% 02.01.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 02.01.2017      Viq output matches zero-yaw output from ViqsNL.m.
%                 Verified C, Dy with complex step using Viq output.
%
% Inputs:
% -------
% Nb              : Number of blades.
% Viz             : Z-component induced velocity (state).
% Vzts            : Local rotorplane flow velocity, zts components.
% rho             : Air density.
% r               : Element radius.
% L               : Element length.
% f               : Prandtl factor.
% Fzts            : Airfoil forces in zts coordinates.
%
% Outputs:
% --------
% Viq             : z,t components of induced velocity.
% C,Dy            : State matrices.

%{
'---Vi_qs---'
Nb
Viz
Vzts
rho
r
L
f
Fzts
%}

C  = zeros(2,1);
Dy = zeros (2,9);

a2 = 1;
CT2 = 1.82;
a1 = 1 - 0.5*sqrt(CT2);
CT1 = 4*a1*(1 - a1);

vfv     = Vzts(1) + f*Viz;
W       = sqrt(vfv^2 + Vzts(2)^2 + Vzts(3)^2);
A       = 2*pi*r*L/Nb;
twopAfW = 2*rho*A*f*W;
val     = 4*pi*rho/(Nb*(twopAfW^2));
rLf     = r*L*f;

Viq = -Fzts(1:2)/twopAfW;

%{
'------'
W
A
Viq
-real((a1/f)*Vzts(1))
%}

if (real(Viq(1)) < -real((a1/f)*Vzts(1)))

   % Employ a Glauert-type formula for induced velocity.
   V = sqrt(Vzts(1)^2 + Vzts(2)^2 + Vzts(3)^2);
   pAV2 = rho*A*(V^2);
   CC = (Fzts(1)/(0.5*(CT2 - CT1)*pAV2) - CT1/(CT2 - CT1))*(a2 - a1) + a1;
   Viq(1) = -Vzts(1)*CC/f;
   val2 = (Fzts(1)*(a2 - a1)/((0.5*(CT2 - CT1)*pAV2)^2))*(CT2 - CT1)*rho*pi/Nb;
   tworL = 2*r*L;

   Dy(1,1)   = -CC/f + (Vzts(1)/f)*val2*tworL*Vzts(1);
   Dy(1,2)   = (Vzts(1)/f)*val2*tworL*Vzts(2);
   Dy(1,3)   = (Vzts(1)/f)*val2*tworL*Vzts(3);
   Dy(1,4)   = (Vzts(1)/f)*val2*L*(V^2);
   Dy(1,5)   = (Vzts(1)/f)*val2*r*(V^2);
   Dy(1,6)   =  CC*Vzts(1)/(f^2);
   Dy(1,7)   = -(Vzts(1)/f)*(a2 - a1)/(0.5*(CT2 - CT1)*pAV2);

   Dy(2,1)   = Fzts(2)*val*rLf*vfv/W;
   Dy(2,2)   = Fzts(2)*val*rLf*Vzts(2)/W;
   Dy(2,3)   = Fzts(2)*val*rLf*Vzts(3)/W;
   Dy(2,4)   = Fzts(2)*val*L*f*W;
   Dy(2,5)   = Fzts(2)*val*r*f*W;
   Dy(2,6)   = Dy(2,1)*Viz + Fzts(2)*val*r*L*W;
   Dy(2,8)   = -1/twopAfW;
   C(2)      = Dy(2,1)*f;

%{
'------'
'Glauert'
V
Ct = Fzts(1)/(0.5*pAV2)
Viq(1)
%}

else

   % Employ the momentum balance formula for induced velocity.
   Dy(:,1)   = Fzts(1:2)*val*rLf*vfv/W;
   Dy(:,2)   = Fzts(1:2)*val*rLf*Vzts(2)/W;
   Dy(:,3)   = Fzts(1:2)*val*rLf*Vzts(3)/W;
   Dy(:,4)   = Fzts(1:2)*val*L*f*W;
   Dy(:,5)   = Fzts(1:2)*val*r*f*W;
   Dy(:,6)   = Dy(:,1)*Viz + Fzts(1:2)*val*r*L*W;
   Dy(:,7:8) = -diag([1;1])/twopAfW;
   C(:,:)    = Dy(:,1)*f;

end

%'---end Viqs---'