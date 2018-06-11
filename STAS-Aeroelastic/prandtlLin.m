function [f,z,Dy] = prandtlLin (Nb,r,Do,Uzts)
%
% Compute the linearization of Prandtl's tip loss/finite blade function.
%
%   States:           y vector:         u vector:
%                     Uzts
%                     r
%                     Do
%                     z
%
% Version:        Changes:
% --------        -------------
% 27.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 27.12.2017      Checked against complex step using f,z output.  Also
%                 checked that f agrees with prandtl.m.
%
% Inputs:
% -------
% Nb              : Number of blades.
% r               : Radius of element.
% Do              : Outer diameter of rotor (deformed blade).
% Uzts            : Local rotorplane windspeed, zts components.
%
% Outputs:
% --------
% f,z             : Prandtl factor and intermediate variable.
% Dy              : State space.

Dy = zeros(2,6);

phi = atan2c (Uzts(1),Uzts(2));
sp  = sin(phi);
tp  = tan(phi);

if (real(sp) < (eps^0.25))

   % Should not be the case in BEM, but nonetheless prevent divide-by-zero.
   z = 0;
   f = 1;
   % Dy = 0.

else

   U2   =  Uzts(1)^2 + Uzts(2)^2;
   dpdz =  Uzts(2)/U2;
   dpdt = -Uzts(1)/U2;

   k = Nb*(r - 0.5*Do)/(2*r);
   kst = k/(sp*tp);

   z = k/sp;
   f = (2/pi)*acos(exp(z));

   Nb2r = Nb/(2*r);

   Dy(2,1) = -kst*dpdz;                      % Uz
   Dy(2,2) = -kst*dpdt;                      % Ut
   Dy(2,5) = -0.5*Nb2r/sp;                   % Do
   Dy(2,4) =  Nb2r*(1 - (r - 0.5*Do)/r)/sp;  % r

   if (real(f) < 0.1)
      f = 0.1;
      % Dy(1,:) = 0;
   else
      Dy(1,6) = -(2/pi)*exp(z)/sqrt(1 - exp(2*z));
   end

end
