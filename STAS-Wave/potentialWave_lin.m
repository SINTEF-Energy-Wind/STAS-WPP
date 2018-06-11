function [p,zeta,vx,ax] = potentialWave_lin (Hs,Tp,d,Dia,rho,f,theta,z,t)
%
% Simulates a wave about a monopile in the time domain.
% Pressures are computed using first-order wave kinematics
% but including the v^2 terms in Bernoulli's equation.
%
% Version:        Changes:
% --------        -------------
% 25.03.2015      Original code.
%
% Version:        Verification:
% --------        -------------
% 25.03.2015      
%
% Inputs:
% -------
% Hs              : Significant wave height.
% Tp              : Characteristic wave period.
% d               : Water depth.
% Dia             : Cylinder diameter.
% rho             : Water density.
% f               : Frequencies to use when discretizing the wave
%                   amplitude spectrum.
% theta           : Location about the cylinder, theta = 0 is
%                   on the lee side relative to the waves.
% z               : Depth coordinates relative to the free surface.
% t               : Vector of times.
%
% Outputs:
% --------
% p               : Ntime*Nth*Nz vector of pressures.
% zeta            : Ntime*Nth vector of surface elevations.
% vx              : Ntime*Nz vector of incident velocities at the
%                   cylinder centerplane.
% ax              : Ntime*Nz vector of incident accelerations at the
%                   cylinder centerplane.

g = 9.807;

Nth = size(theta,1);
Nz = size(z,1);
Nf = size(f,1);
Ntime = size(t,1);

zeta = zeros(Ntime*Nth,1);
p = zeros(Ntime*Nth*Nz,1);
vx = zeros(Ntime*Nz,1);
ax = zeros(Ntime*Nz,1);

% Dispersion relationship and decomposition of wave spectra
% into discrete amplitudes.
[zeta0,eps,k,df] = waveFourier (d,Hs,Tp,f);

% Complex amplitude of velocity for each frequency bin.
[phi0,vt0,vz0,vx0] = waveVelocity (f,zeta0,k,eps,Dia,d,theta,z);

% Matrix of cyclic phases.  Want f row-wise.
eiwt = exp(i*2*pi*t*f');

% First get the wave elevation above the theta stations.
for ith = 1:Nth
   indt = Ntime*(ith-1);
   zeta(indt+1:indt+Ntime) = ...
         real(eiwt*(exp(i*(-0.5*k*Dia*cos(theta(ith)) + eps)).*zeta0));
end

fid = fopen('zeta.txt','w');
for itime = 1:Ntime
   fprintf(fid,'%+5.4e %+5.4e %+5.4e %+5.4e\n',      ...
           t(itime),zeta(itime),zeta(9*Ntime+itime), ...
           zeta(18*Ntime+itime));
end
fclose(fid);

for iz = 1:Nz
   indzf = Nf*Nth*(iz-1);
   indzt = Ntime*Nth*(iz-1);
   indx0 = Nf*(iz-1);
   indx = Ntime*(iz-1);

   for ith = 1:Nth
      indf = indzf + Nf*(ith-1);
      indt = indzt + Ntime*(ith-1);
      izet = Ntime*(ith-1);

      % Get arrays at the current point and all times.
      dphidt = real(i*2*pi*eiwt*(f.*phi0(indf+1:indf+Nf)));
      vt = real(eiwt*vt0(indf+1:indf+Nf));
      vz = real(eiwt*vz0(indf+1:indf+Nf));

      % Compute the pressure at the current point and all times.
      % In the COMPUTATIONAL domain, the hydrostatic pressure
      % does not change.  The dynamic free-surface boundary
      % condition ensures that the surface fluctuations in gz
      % pressure are captured.
      p(indt+1:indt+Ntime) = -rho*dphidt; % -rho*(dphidt + 0.5*(vt.^2 + vz.^2));
%if (iz == Nz)
%[z(iz) theta(ith) p(indt+129)]
%end
   end

   % Compute the incident x-direction flow velocity.
   vx(indx+1:indx+Ntime) = real(eiwt*vx0(indx0+1:indx0+Nf));
   ax(indx+1:indx+Ntime) = real(eiwt*(i*2*pi*f.*vx0(indx0+1:indx0+Nf)));

end

