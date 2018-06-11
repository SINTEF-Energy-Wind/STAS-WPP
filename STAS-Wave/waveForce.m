function F = waveForce (Hs,Tp,Vc,d,rho,f,theta,z,t,znod,Dnod)
%
% This computes the net X-direction nodal forces on a monopile
% due to ocean waves.
%
% Version:        Changes:
% --------        -------------
% 25.03.2015      Original code.
%
% Version:        Verification:
% --------        -------------
% 25.03.2015      Memo AN 15.12.19
%
% Inputs:
% -------
% Hs              : Significant wave height.
% Tp              : Characteristic wave period.
% Vc              : X (wave direction) and Y (transverse) current
%                   velocities in m/s, at each computational
%                   element.
% d               : Water depth.
% rho             : Water density.
% f               : Frequencies to use when discretizing the wave
%                   amplitude spectrum.
% theta           : Location about the cylinder, theta = 0 is
%                   on the lee side relative to the waves.
% z               : Depth coordinates relative to the free surface.
% t               : Vector of times.
% znod            : elevations of the nodes relative to the mean 
%                   free surface, negative = submerged.
% Dnod            : Cylinder diameter at each node.
%
% Outputs:
% --------
% F               : Time history of wave-direction forces on each
%                   node.  Vector of length 2*Ntime*Nnod.  The first
%                   Ntime*Nnod entries are the X forces, and the
%                   second are the Y forces.

Nth = size(theta,1);
Nz = size(z,1);
Nf = size(f,1);
Ntime = size(t,1);
Nnod = size(znod,1);

% Assumes symmetry.
dth = zeros(Nth,1);
dth(1) = 0.5*(theta(2) - theta(1));
dth(Nth) = 0.5*(theta(Nth) - theta(Nth-1));
dth(2:Nth-1) = 0.5*(theta(3:Nth) - theta(1:Nth-2));

% No assumption of symmetry.
%dth = zeros(Nth,1);
%dth(1) = 0.5*(theta(2) - (theta(Nth) - 2*pi));
%dth(Nth) = 0.5*((theta(1) + 2*pi) - theta(Nth-1));
%dth(2:Nth-1) = 0.5*(theta(3:Nth) - theta(1:Nth-2));

dz = zeros(Nz,1);
dz(1) = z(2) - z(1);
dz(Nz) = z(Nz) - z(Nz-1);
dz(2:Nz-1) = 0.5*(z(3:Nz) - z(1:Nz-2));

% The fraction of the element length at the node.
fz = zeros(Nz,1);
fz(1) = 0.5;
fz(Nz) = 0.5;
fz(2:Nz-1) = 0.5*(z(2:Nz-1) - z(1:Nz-2))./dz(2:Nz-1);

F = sparse(2*Ntime*Nnod,1);

% Get the time histories of fluctuations in pressure and elevation.
% Assume a reference diameter here, the pressure will later act
% over the appropriate diameter.  OK for small variations in 
% diameter.
Dref = mean(Dnod);
[p,zeta,vx,ax] = potentialWave (Hs,Tp,d,Dref,rho,f,theta,z,t);

% For each unstretched element on the (computational) cylinder
% surface ...
for iz = 1:Nz
   indzt = Ntime*Nth*(iz-1);
   indx = Ntime*(iz-1);

   for ith = 1:Nth
      indt = indzt + Ntime*(ith-1);
      izet = Ntime*(ith-1);
      
      % Stretch the element to an effective z coordinate on the
      % real cylinder.
      zp = ((z(iz) + d)*(d + zeta(izet+1:izet+Ntime))/d) - d;
      dzp = dz(iz)*(d + zeta(izet+1:izet+Ntime))/d;
      Diap = interp1(znod,Dnod,zp,'extrap');

      % dF = -p dA' cos(theta).  Force per unit length.
      dF = -p(indt+1:indt+Ntime).*0.5.*Diap*dth(ith)*cos(theta(ith));

      % Distribute to the appropriate nodes.
      for ii = 1:Nnod
         ind = Ntime*(ii-1);

         C = zeros(Ntime,1);

         if (ii ~= Nnod)
            a = -1/(znod(ii+1) - znod(ii));
            b = znod(ii+1)/(znod(ii+1) - znod(ii));

            LB = min(max(zp-fz(iz)*dzp,znod(ii)),znod(ii+1));
            RB = min(max(zp+(1-fz(iz))*dzp,znod(ii)),znod(ii+1));

            C = C + 0.5*a*(RB.^2 - LB.^2) + b*(RB - LB);
 
         end

         if (ii ~= 1)
            a = 1/(znod(ii) - znod(ii-1));
            b = -znod(ii-1)/(znod(ii) - znod(ii-1));

            LB = min(max(zp-fz(iz)*dzp,znod(ii-1)),znod(ii));
            RB = min(max(zp+(1-fz(iz))*dzp,znod(ii-1)),znod(ii));

            C = C + 0.5*a*(RB.^2 - LB.^2) + b*(RB - LB);
 
         end

         % Symmetry: double the force increment.
% (Comment out for Morison.)
         F(ind+1:ind+Ntime) = F(ind+1:ind+Ntime) + 2*C.*dF;

         % No symmetry.
         %F(ind+1:ind+Ntime) = F(ind+1:ind+Ntime) + C.*dF;

      end

   end

   % Compute the viscous drag on the stretched z location.
   % Stretch the profile.  Cd = 1 for design.
   % dFd is the force per unit length.
   Cd = 1;
%Cd = 0.6;  % For comparing against laboratory data.
   Velx = vx(indx+1:indx+Ntime) + Vc(1,iz);
   Vmag = sqrt(Velx.^2 + Vc(2,iz)^2);
   qi = ceil(interp1(theta,[1:Nth],pi/2));
   iztc = Ntime*(qi-1);
   zp = ((z(iz) + d)*(d + zeta(iztc+1:iztc+Ntime))/d) - d;
   dzp = dz(iz)*(d + zeta(izet+1:izet+Ntime))/d;
   Diap = interp1(znod,Dnod,zp,'extrap');
   dFdx = Cd*0.5*rho*Diap.*Vmag.*Velx;  % (Comment out for Morison.)
   dFdy = Cd*0.5*rho*Diap.*Vmag.*Vc(2);

% Uncomment for Morison.  Also see lines above to comment out.
%dFdx = Cd*0.5*rho*Diap.*Vmag.*Velx + 2*rho*(0.25*pi*(Diap.^2)).*ax(indx+1:indx+Ntime);

   % Distribute to the appropriate nodes.
   ioff = Ntime*Nnod;
   for ii = 1:Nnod
      ind = Ntime*(ii-1);

      C = zeros(Ntime,1);

      if (ii ~= Nnod)
         a = -1/(znod(ii+1) - znod(ii));
         b = znod(ii+1)/(znod(ii+1) - znod(ii));

         LB = min(max(zp-fz(iz)*dzp,znod(ii)),znod(ii+1));
         RB = min(max(zp+(1-fz(iz))*dzp,znod(ii)),znod(ii+1));

         C = C + 0.5*a*(RB.^2 - LB.^2) + b*(RB - LB);

      end

      if (ii ~= 1)
         a = 1/(znod(ii) - znod(ii-1));
         b = -znod(ii-1)/(znod(ii) - znod(ii-1));

         LB = min(max(zp-fz(iz)*dzp,znod(ii-1)),znod(ii));
         RB = min(max(zp+(1-fz(iz))*dzp,znod(ii-1)),znod(ii));

         C = C + 0.5*a*(RB.^2 - LB.^2) + b*(RB - LB);

      end

      F(ind+1:ind+Ntime) = F(ind+1:ind+Ntime) + C.*dFdx;

      F(ioff+ind+1:ioff+ind+Ntime) = F(ioff+ind+1:ioff+ind+Ntime) ...
                                   + C.*dFdy;

   end

end

