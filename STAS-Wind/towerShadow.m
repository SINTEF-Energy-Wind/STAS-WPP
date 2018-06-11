function [SzzR,SzzI] = towerShadow (Neb,r,x,Dt,df,W,Vzrp,N)
%
% Adds spikes to the Szz spectrum in order to account for the effects of
% the tower on excitation of blade vibration.  Potential theory is used,
% employing Burton p 234, Eq 5.21.
%
% Presently only the uz (axial) component of windspeed is modelled.
%
% Version:        Changes:
% --------        -------------
% 30.08.2015      Adapted from the towerDam subroutine in the
%                 modalAnalysis.f90 Fortran code.
% 03.08.2017      Adapted for complex step derivatives.
%
% Version:        Verification:
% --------        -------------
% 30.08.2015      
% 03.08.2017      
%
% Inputs:
% -------
% Neb             : Number of elements per blade.
% r               : Radial coordinate of each element (size Neb).
% x               : Distance from tower centerline to element (size Neb).
% Dt              : Tower diameter at height of element (size Neb).
% df              : Frequency bin width.
% W               : Rotor rotational speed.
% Vzrp            : Nominal windspeed at the rotorplane.  May include 
%                   induction.
% N               : Number of Fourier terms to use.
%
% Outputs:
% --------
% Szz             : A matrix of velocity auto- and cross-spectra, size
%                   N-by-(3*Neb)^2.

SzzR = sparse(N,(3*Neb)^2);
SzzI = sparse(N,(3*Neb)^2);

Nsteps = 3600;

dpsi = 2*pi/Nsteps;
psi = dpsi*[0:(Nsteps-1)].';
cp = cos(psi + pi/2);  % pi/2 so that 0 deg aligns with the X^r axis.
sp = sin(psi + pi/2);

% Mean offset.
c0 = zeros(Neb,1);
vec = ones(Nsteps,1);
uz = zeros(Nsteps,Neb);
for iel = 1:Neb
   y = sqrt((r(iel)*sp).^2 + (r(iel)*(1 - cp)).^2);
   uz(:,iel) = -Vzrp(iel)*((0.5d0*Dt(iel))^2)*(x(iel)^2 - y.^2) ...
     ./ ((x(iel)^2 + y.^2).^2);
   c0(iel) = (uz(:,iel).')*vec/Nsteps;
end

cR = zeros(N,3*Neb);
cI = zeros(N,3*Neb);
e2p3 = exp(i*[1:N]'*2*pi/3);
e4p3 = exp(i*[1:N]'*4*pi/3);
for ifreq = 1:N

   e = exp(-i*ifreq*psi);

   % Generate the coefficients for one blade.
   for iel = 1:Neb
      cR(ifreq,iel) = (uz(:,iel).')*real(e)/Nsteps;
      cI(ifreq,iel) = (uz(:,iel).')*imag(e)/Nsteps;
   end

   % Coefficients for other blades, apply an offset.
   [cR(ifreq,Neb+1:2*Neb),cI(ifreq,Neb+1:2*Neb)] =     ...
         cmult(cR(ifreq,1:Neb),cI(ifreq,1:Neb),        ...
               real(e2p3(ifreq)),imag(e2p3(ifreq)));
   [cR(ifreq,2*Neb+1:3*Neb),cI(ifreq,2*Neb+1:3*Neb)] = ...
         cmult(cR(ifreq,1:Neb),cI(ifreq,1:Neb),        ...
               real(e4p3(ifreq)),imag(e4p3(ifreq)));

end

% The complex spectra.
for iel1 = 1:3*Neb
   for iel2 = 1:3*Neb
      icols = 3*Neb*(iel1-1) + iel2;

      [SzzR(:,icols),SzzI(:,icols)] =       ... 
            cdotmult(cR(:,iel1),cI(:,iel1), ...
                     cR(:,iel2),-cI(:,iel2));

   end
end

SzzR = SzzR/df;
SzzI = SzzI/df;
