function [phi,vt,vz,vx] = waveVelocity (f,zeta0,k,eps,Dia,d,theta,z)
%
% This function computes the complex amplitudes of the potential,
% and the tangential and vertical velocity components at the
% surface of a cylinder in regular waves.
%
% Version:        Changes:
% --------        -------------
% 24.03.2015      Original code.
%
% Version:        Verification:
% --------        -------------
% 24.03.2015      
%
% Inputs:
% -------
% f               : Frequencies.
% zeta0           : Wave amplitudes, one per frequency.
% k               : Wave number, one per frequency.
% eps             : Random phases, one per frequency.
% Dia             : Cylinder diameter.  Assumed to be a scalar.
% d               : Water depth.
% theta           : Location about the cylinder, theta = 0 is
%                   on the lee side relative to the waves.
% z               : Depth coordinates relative to the free surface.
% Nm              : Number of series terms to retain.
%
% Outputs:
% --------
% phi             : Complex amplitude of the velocity potential.
% vt,vz           : Complex amplitudes of the theta- and z-
%                   direction velocities at the given points on
%                   the cylinder surface.
% vx              : Complex amplitude of the INCIDENT x-direction
%                   velocity at the centerplane of the cylinder.
%                   

g = 9.807;
R = 0.5*Dia;
omega = 2*pi*f;
kR = k*R;
kd = k*d;

Nf = size(f,1);
Nz = size(z,1);
Nt = size(theta,1);
Nm = 10;

m = [0:Nm-1]';

phi = zeros(Nf*Nt*Nz,1);
vt  = zeros(Nf*Nt*Nz,1);
vz  = zeros(Nf*Nt*Nz,1);
vx  = zeros(Nf*Nz,1);

% Terms which are not dependent on z, t, or m, only f.
gzw = g*zeta0./omega;
eieps = exp(i*eps);

for im = 1:Nm

   if (m(im) ~= 0)
      em = 2;
   else
      em = 1;
   end

   % Depend only on m, f.
   iepst = (i^(m(im)+1))*((-1)^(m(im)+1))*em*m(im)*sin(m(im)*theta);
   iepsv = (i^(m(im)+1))*((-1)^m(im))*em*cos(m(im)*theta);

    jmR = besselj(m(im),  kR);
   jm1R = besselj(m(im)+1,kR);
    ymR = bessely(m(im),  kR);
   ym1R = bessely(m(im)+1,kR);

   jj = -k.*jm1R + (m(im)/R)*jmR;
   yy = -k.*ym1R + (m(im)/R)*ymR;
   fmR = jj./(jj - i*yy);

   for iz = 1:Nz
      indz = (Nf*Nt)*(iz-1);

      kzd = k*(z(iz) + d);
      csh = (exp(kzd - kd) + exp(-(kzd + kd)))./(1 + exp(-2*kd));
       sh = (exp(kzd - kd) - exp(-(kzd + kd)))./(1 + exp(-2*kd));

      jfy = (jmR - fmR.*(jmR - i*ymR));
       ge = gzw.*eieps;
       pp = ge.*csh.*jfy;
      vvt = pp/R;
      vvz = k.*ge.*sh.*jfy;

      for it = 1:Nt
         ind = indz + Nf*(it-1);

         phi(ind+1:ind+Nf) = phi(ind+1:ind+Nf) +  pp*iepsv(it);
          vt(ind+1:ind+Nf) =  vt(ind+1:ind+Nf) + vvt*iepst(it);
          vz(ind+1:ind+Nf) =  vz(ind+1:ind+Nf) + vvz*iepsv(it);
%[theta(it) z(iz) (vvz*iepsv(it))(60)]
%if (iz == 6) && (it == 1)
%[m(im) phi(ind+60)]
%end
      end

   end

end

% Compute the undisturbed velocity on the z section, at a
% location corresponding to the cylinder centerplane, x=0.
for iz = 1:Nz
   indz = Nf*(iz-1);
   kzd = k*(z(iz) + d);
   csh = (exp(kzd - kd) + exp(-(kzd + kd)))./(1 + exp(-2*kd));
   vx(indz+1:indz+Nf) = k.*gzw.*csh.*eieps;
end


