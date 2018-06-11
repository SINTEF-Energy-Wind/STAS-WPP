function [SpsiR,SpsiI] = towerShadowMBC (Neb,r,x,Dt,df,W,Vzrp,N)
%
% Adds spikes to the Spsi spectrum in order to account for the effects of
% the tower on excitation of blade vibration.  Potential theory is used,
% employing Burton p 234, Eq 5.21.
%
% Presently only the uz (axial) component of windspeed is modelled.  It
% is converted to multi-blade coordinates before obtaining the spectra
% numerically.
%
% Version:        Changes:
% --------        -------------
% 01.09.2015      Original code.
% 03.08.2017      Adapted for complex step derivatives.
%
% Version:        Verification:
% --------        -------------
% 01.09.2015      
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
% Spsi            : A matrix of velocity auto- and cross-spectra, size
%                   N-by-(3*Neb)^2.

SpsiR = sparse(N,(3*Neb)^2);
SpsiI = sparse(N,(3*Neb)^2);

Nsteps = 3600;

% The formula for tower shadow is expressed with zero angle defined
% as blade aligned with tower.  The coordinate system I am using has
% psi = 0 defined as parallel with the X^r axis, which is pi/2 offset
% from the tower.  Thus for terms used in the formulas for uz, add
% pi/2 to the cosine and sine.
dpsi = 2*pi/Nsteps;
psi = dpsi*[0:(Nsteps-1)].';
cp = cos(psi);
sp = sin(psi);
cp2 = cp*cos(2*pi/3) - sp*sin(2*pi/3);
cp4 = cp*cos(4*pi/3) - sp*sin(4*pi/3);
sp2 = sp*cos(2*pi/3) + cp*sin(2*pi/3);
sp4 = sp*cos(4*pi/3) + cp*sin(4*pi/3);

ct  = cp *cos(pi/2) - sp *sin(pi/2);
ct2 = cp2*cos(pi/2) - sp2*sin(pi/2);
ct4 = cp4*cos(pi/2) - sp4*sin(pi/2);
st  = sp *cos(pi/2) + cp *sin(pi/2);
st2 = sp2*cos(pi/2) + cp2*sin(pi/2);
st4 = sp4*cos(pi/2) + cp4*sin(pi/2);

% Mean offset.
c0 = zeros(3*Neb,1);
vec = ones(Nsteps,1);
uz = zeros(Nsteps,3);
upsi = zeros(Nsteps,3*Neb);
y    = zeros(Nsteps,3);
for iel = 1:Neb
   y(:,1) = sqrt((r(iel)*st ).^2 + (r(iel)*(1 - ct )).^2);
   y(:,2) = sqrt((r(iel)*st2).^2 + (r(iel)*(1 - ct2)).^2);
   y(:,3) = sqrt((r(iel)*st4).^2 + (r(iel)*(1 - ct4)).^2);
   uz(:,1) = -Vzrp(iel)*((0.5*Dt(iel))^2)*(x(iel)^2 - y(:,1).^2) ...
          ./ ((x(iel)^2 + y(:,1).^2).^2);
   uz(:,2) = -Vzrp(iel)*((0.5*Dt(iel))^2)*(x(iel)^2 - y(:,2).^2) ...
          ./ ((x(iel)^2 + y(:,2).^2).^2);
   uz(:,3) = -Vzrp(iel)*((0.5*Dt(iel))^2)*(x(iel)^2 - y(:,3).^2) ...
          ./ ((x(iel)^2 + y(:,3).^2).^2);
   upsi(:,iel)       =   (uz(:,1)     + uz(:,2)      + uz(:,3)     )/3;
   upsi(:,iel+Neb)   = 2*(uz(:,1).*cp + uz(:,2).*cp2 + uz(:,3).*cp4)/3;
   upsi(:,iel+2*Neb) = 2*(uz(:,1).*sp + uz(:,2).*sp2 + uz(:,3).*sp4)/3;
   c0(iel)       = upsi(:,iel).'      *vec/Nsteps;
   c0(iel+Neb)   = upsi(:,iel+Neb).'  *vec/Nsteps;
   c0(iel+2*Neb) = upsi(:,iel+2*Neb).'*vec/Nsteps;

end

cR = zeros(N,3*Neb);
cI = zeros(N,3*Neb);
for ifreq = 1:N

%   e = exp(-i*ifreq*psi);
   e = exp(-i*3*ifreq*psi);

   % Generate the coefficients.
   for iel = 1:Neb

      cR(ifreq,iel)       = upsi(:,iel).'      *real(e)/Nsteps;
      cR(ifreq,iel+Neb)   = upsi(:,iel+Neb).'  *real(e)/Nsteps;
      cR(ifreq,iel+2*Neb) = upsi(:,iel+2*Neb).'*real(e)/Nsteps;

      cI(ifreq,iel)       = upsi(:,iel).'      *imag(e)/Nsteps;
      cI(ifreq,iel+Neb)   = upsi(:,iel+Neb).'  *imag(e)/Nsteps;
      cI(ifreq,iel+2*Neb) = upsi(:,iel+2*Neb).'*imag(e)/Nsteps;

   end

end

% The complex spectra.
for iel1 = 1:3*Neb
   for iel2 = 1:3*Neb
      icols = 3*Neb*(iel1-1) + iel2;

      [SpsiR(:,icols),SpsiI(:,icols)] =    ...
         cdotmult(cR(:,iel1),cI(:,iel1),   ...
                  cR(:,iel2),-cI(:,iel2));

   end
end

SpsiR = SpsiR/df;
SpsiI = SpsiI/df;
