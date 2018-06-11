function [SpsiR,SpsiI] = windShearMBC (Neb,r,H,z0,df,W,Vinf,N)
%
% "Standard" wind shear is a logarithmic profile of mean windspeed with
% elevation.  When a point is swept around this mean windspeed profile,
% a Fourier decomposition can be used to represent the perturbation in
% the observed relative windspeed.  First, the observed windspeed is 
% converted to multi-blade coordinates.  The coefficients on the Fourier
% decomposition are then computed numerically.  The resulting (complex)
% cross-spectra are output in the Spsi matrix.
%
% Version:        Changes:
% --------        -------------
% 01.09.2015      Original code.
% 02.08.2017      Adapted for complex step derivatives.
%
% Version:        Verification:
% --------        -------------
% 01.09.2015      
% 02.08.2017      Derivatives verified against finite difference.
%                 Reproduces the values from the previous version.
%
% Inputs:
% -------
% Neb             : Number of elements per blade.
% r               : Radial coordinate of each element.
% H               : Hub height.
% z0              : Surface roughness length.
% df              : Frequency bin width.
% W               : Rotor rotational speed.
% Vinf            : Mean windspeed at hub height.
% N               : Number of Fourier terms to use.
%
% Outputs:
% --------
% Spsi            : A matrix of velocity auto- and cross-spectra, size
%                   N-by-(3*Neb)^2.

SpsiR = sparse(N,(3*Neb)^2);
SpsiI = sparse(N,(3*Neb)^2);

Nsteps = 72;

dpsi = 2*pi/Nsteps;
psi = dpsi*[0:(Nsteps-1)].';
cp  = cos(psi);
sp  = sin(psi);
cp2 = cp*cos(2*pi/3) - sp*sin(2*pi/3);
cp4 = cp*cos(4*pi/3) - sp*sin(4*pi/3);
sp2 = sp*cos(2*pi/3) + cp*sin(2*pi/3);
sp4 = sp*cos(4*pi/3) + cp*sin(4*pi/3);

% Mean offset.
c0   = zeros(3*Neb,1);
vec  =  ones(Nsteps,1);
uz   = zeros(Nsteps,3);
upsi = zeros(Nsteps,3*Neb);
y    = zeros(Nsteps,3);
for iel = 1:Neb
   y(:,1) = H + r(iel)*sp;  % psi = 0 is the X^r axis.
   y(:,2) = H + r(iel)*sp2;
   y(:,3) = H + r(iel)*sp4;
   uz(:,1) = Vinf*(log(y(:,1)/z0)/log(H/z0) - 1);
   uz(:,2) = Vinf*(log(y(:,2)/z0)/log(H/z0) - 1);
   uz(:,3) = Vinf*(log(y(:,3)/z0)/log(H/z0) - 1);
   upsi(:,iel)       =   (uz(:,1)     + uz(:,2)      + uz(:,3)     )/3;
   upsi(:,iel+Neb)   = 2*(uz(:,1).*cp + uz(:,2).*cp2 + uz(:,3).*cp4)/3;
   upsi(:,iel+2*Neb) = 2*(uz(:,1).*sp + uz(:,2).*sp2 + uz(:,3).*sp4)/3;
   c0(iel)       = upsi(:,iel).'      *vec/Nsteps;
   c0(iel+Neb)   = upsi(:,iel+Neb).'  *vec/Nsteps;
   c0(iel+2*Neb) = upsi(:,iel+2*Neb).'*vec/Nsteps;

end

cR = zeros(N,3*Neb);
cI = zeros(N,3*Neb);
for ifreq = 1:N  % Blade passing frequency.

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

% Should iel1 or iel2 be conjugate?
      [SpsiR(:,icols),SpsiI(:,icols)] =    ...
         cdotmult(cR(:,iel1),cI(:,iel1),   ...
                  cR(:,iel2),-cI(:,iel2));
%      [SpsiR(:,icols),SpsiI(:,icols)] =    ...
%         cdotmult(cR(:,iel1),cI(:,iel1),   ...
%                  cR(:,iel2),-cI(:,iel2));

   end
end

SpsiR = SpsiR/df;
SpsiI = SpsiI/df;

