function [Viy,C,Dy] = Viyaw (Viz,Viq,Vzts,r,Do,f,W)
%
% Quasi-steady induced velocity in yaw.
%
%   States:           y vector:         u vector:
%   Viz               Vzts   1:3
%                     Viq    4,5
%                     r       6
%                     Do      7
%                     f       8
%                     W       9
%
% Version:        Changes:
% --------        -------------
% 02.01.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 02.01.2017      Viy matches output from ViqsNL.m.  Checked using complex
%                 step with Viy output.
%
% Inputs:
% -------
% Viz             : z-component induced velocity (state).
% Viq             : z,t quasi-steady induced velocity.
% Vzts            : Local rotorplane flow velocity, zts components.
% r               : Element radius.
% Do              : Outer diameter (deformed blade).
% f               : Prandtl factor.
% W               : Local flow velocity magnitude.
%
% Outputs:
% --------
% Viy             : z,t components of induced velocity in yaw.
% C,Dy            : State matrices.

C  = zeros (2,1);
Dy = zeros (2,9);

vfv     = Vzts(1) + f*Viz;
rD      = r/Do;
mag     = sqrt(Vzts(2)^2 + Vzts(3)^2);

if (real(mag) > (sqrt(eps)*abs(real(Vzts(1)))))
   z     = vfv/W;
   tcz   = tan(0.5*acos(z));
   zsz   = 1/((z + 1)*sqrt(1 - z^2));
   kap   = 1 + 2*rD*tcz*Vzts(3)/mag;
   dkdvz = -2*rD*Vzts(3)*zsz/(W*mag);
   dkdvt = -2*rD*tcz*Vzts(2)*Vzts(3)/(mag^3);
   dkdvs = 2*rD*tcz*(Vzts(2)^2)/(mag^3);
   dkdvi = f*dkdvz;
   dkdf  = Viz*dkdvz;
   dkdW  = 2*rD*Vzts(3)*zsz*vfv/((W^2)*mag);
   dkdr  = (2/Do)*tcz*Vzts(3)/mag;
   dkdD  = -(2*r/(Do^2))*tcz*Vzts(3)/mag;
   Viy   = kap*Viq;
   Dy(:,1) = Viq*dkdvz;
   Dy(:,2) = Viq*dkdvt;
   Dy(:,3) = Viq*dkdvs;
   Dy(1,4) = kap;
   Dy(2,5) = kap;
   Dy(:,6) = Viq*dkdr;
   Dy(:,7) = Viq*dkdD;
   Dy(:,8) = Viq*dkdf;
   Dy(:,9) = Viq*dkdW;   
   C(:,:)  = Viq*dkdvi;
else  % Infinitesimal yaw angle.
   Viy = Viq;
   Dy(:,4:5) = eye(2);
end

