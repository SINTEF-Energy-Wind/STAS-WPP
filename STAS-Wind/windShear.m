function [SzzR,SzzI] = windShear (Neb,r,H,z0,df,W,Vinf,N)
%
% "Standard" wind shear is a logarithmic profile of mean windspeed with
% elevation.  When a point is swept around this mean windspeed profile,
% a Fourier decomposition can be used to represent the perturbation in
% the observed relative windspeed.  The coefficients on the Fourier
% decomposition are here computed numerically.  The resulting (complex)
% cross-spectra are output in the Szz matrix.
%
% Version:        Changes:
% --------        -------------
% 25.08.2015      Adapted from the windShear subroutine in the
%                 modalAnalysis.f90 Fortran code.
% 02.08.2017      Adapted for complex step derivatives.  This involves
%                 separating the real and imaginary components of Szz,
%                 such that EACH of these can be complex, where the
%                 complex part is now associated with the derivative.
%
% Version:        Verification:
% --------        -------------
% 25.08.2015      
% 02.08.2017      Matches the windShear.m from STAS Wind 0.1.  Complex
%                 step verified against finite difference.
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
% SzzR,SzzI       : Real and imaginary components of a matrix of
%                   velocity auto- and cross-spectra, size N-by-
%                   (3*Neb)^2.

SzzR = zeros(N,(3*Neb)^2);
SzzI = zeros(N,(3*Neb)^2);

Nsteps = 72;

dpsi = 2*pi/Nsteps;
psi = dpsi*[0:(Nsteps-1)].';
sp = sin(psi);

% Mean offset.
c0 = zeros(Neb,1);
vec = ones(Nsteps,1);
uz = zeros(Nsteps,Neb);
for iel = 1:Neb
   y = H + r(iel)*sp;  % psi = 0 is the X^r axis.
   uz(:,iel) = Vinf*(log(y/z0)/log(H/z0) - 1);
   c0(iel) = (uz(:,iel).')*vec/Nsteps;
end

cR = zeros(N,3*Neb);
cI = zeros(N,3*Neb);
e2p3 = exp(i*[1:N]'*2*pi/3);
e4p3 = exp(i*[1:N]'*4*pi/3);
for ifreq = 1:N

   e = exp(-i*ifreq*psi);  % None of these can be perturbed by
                           % the inputs.

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

% Where does the conjugate belong, on iel1 or iel2?
      [SzzR(:,icols),SzzI(:,icols)] =       ... 
            cdotmult(cR(:,iel1),cI(:,iel1), ...
                     cR(:,iel2),-cI(:,iel2));
%      [SzzR(:,icols),SzzI(:,icols)] =        ... 
%            cdotmult(cR(:,iel1),-cI(:,iel1), ...
%                     cR(:,iel2),cI(:,iel2));

   end
end

SzzR = SzzR/df;
SzzI = SzzI/df;
