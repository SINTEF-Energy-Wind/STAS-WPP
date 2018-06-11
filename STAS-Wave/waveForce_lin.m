function F = waveForce_lin (Hs,Tp,d,Dia,rho,f,theta,z,t,znod)
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
% znod            : elevations of the nodes relative to the mean 
%                   free surface, negative = submerged.
%
% Outputs:
% --------
% F               : Time history of wave-direction forces on each
%                   node.  Vector of length Ntime*Nnod.

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

F = sparse(Ntime*Nnod,1);

% Get the time histories of fluctuations in pressure and elevation.
[p,zeta,vx,ax] = potentialWave_lin (Hs,Tp,d,Dia,rho,f,theta,z,t);

% For each unstretched element on the (computational) cylinder
% surface ...


fid5 = fopen('out.txt','w');

for iz = 1:Nz
   indzt = Ntime*Nth*(iz-1);
   indx = Ntime*(iz-1);

fp1 = zeros(Ntime,1);
%fp2 = zeros(Ntime,1);
%if(iz == 20)
%fid4 = fopen('fth.txt','w');
%end

   for ith = 1:Nth
      indt = indzt + Ntime*(ith-1);
      izet = Ntime*(ith-1);
      
      % Stretch the element to an effective z coordinate on the
      % real cylinder.
      zp = z(iz)*ones(Ntime,1); % ((z(iz) + d)*(d + zeta(izet+1:izet+Ntime))/d) - d;
      dzp = dz(iz)*ones(Ntime,1); % dz(iz)*(d + zeta(izet+1:izet+Ntime))/d;

      % dF = -p dA' cos(theta).  Force per unit length.
      dF = -p(indt+1:indt+Ntime).*0.5*Dia*dth(ith)*cos(theta(ith));

fp1 = fp1 + dF;
%fp2 = fp2 - dz(iz)*p(indt+1:indt+Ntime)*0.5*Dia*dth(ith)*cos(theta(ith));

%if (iz == 20)
%[theta(ith) cos(theta(ith))]
%[zeta(izet+1:izet+20) (d + zeta(izet+1:izet+20))/d ...
%dF(1:20) (-dz(iz)*p(indt+1:indt+Ntime)*0.5*Dia*dth(ith)*cos(theta(ith)))(1:20)]
%end
%if (iz == 20)
%fprintf(fid4,'%+5.4e %+5.4e %+5.4e %+5.4e %+5.4e %+5.4e\n', ...
%zeta(izet+9),(d + zeta(izet+9))/d, ...
%dF(9),fp1(9), ...
%(-dz(iz)*p(indt+1:indt+Ntime)*0.5*Dia*dth(ith)*cos(theta(ith)))(9),fp2(9));
%end

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
         F(ind+1:ind+Ntime) = F(ind+1:ind+Ntime) + 2*C.*dF;

         % No symmetry.
         %F(ind+1:ind+Ntime) = F(ind+1:ind+Ntime) + C.*dF;

if (ii == 4) % && (iz == Nz-1)
%fprintf(fid5,'%+5.4e %+5.4e %+5.4e %+5.4e %+5.4e %+5.4e %+5.4e %+5.4e\n', ...
%        z(iz),zp(129),dzp(129),theta(ith),C(129),dF(129),C(129)*dF(129),F(ind+129));
fprintf(fid5,'%+5.4e %+5.4e %+5.4e %+5.4e %+5.4e %+5.4e %+5.4e %+5.4e\n', ...
        z(iz),zp(132),dzp(132),theta(ith),C(132),dF(132),C(132)*dF(132),F(ind+132));
%fprintf(fid5,'%+5.4e %+5.4e %+5.4e %+5.4e %+5.4e %+5.4e %+5.4e %+5.4e\n', ...
%        z(iz),zp(74),dzp(74),theta(ith),C(74),dF(74),C(74)*dF(74),F(ind+74));

end

      end


%[zeta(izet+1:izet+200) qind(1:200) F(Ntime*20 + [1:200]') F(Ntime*19 + [1:200]')]

%if (iz == Nz) && (ith == 1)
%[t zp qind qup qdn]
%end

   end

%if (iz == 20)
%fclose(fid4);
%end

   % Compute the viscous drag on the stretched z location.
   % Stretch the profile.  Cd = 1.
   % dFd is the force per unit length.
%   Cd = 1;
Cd = 0.6;
   qi = ceil(interp1(theta,[1:Nth],pi/2));
   iztc = Ntime*(qi-1);
   zp = ((z(iz) + d)*(d + zeta(iztc+1:iztc+Ntime))/d) - d;
   dzp = dz(iz)*(d + zeta(izet+1:izet+Ntime))/d;
   dFd = Cd*0.5*rho*Dia.*(vx(indx+1:indx+Ntime).^2).*sign(vx(indx+1:indx+Ntime));


% Full Morison (inertia only).
%dFd = 0.25*2*pi*rho*(Dia^2)*ax(indx+1:indx+Ntime);




dFd(:) = 0;

if (iz == 60)
izet = Ntime*9;
fid3 = fopen('F1.txt','w');
for it = 1:Nf
   fprintf(fid3,'%+5.4e %+5.4e %+5.4e %+5.4e %+5.4e %+5.4e\n', ...
           t(it),zeta(izet+it),(d + zeta(izet+it))/d,vx(indx+it),fp1(it),dFd(it));
end
fclose(fid3);
%[t p(indt+1:indt+Ntime)]
end



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

      F(ind+1:ind+Ntime) = F(ind+1:ind+Ntime) + C.*dFd;

   end

end

fclose(fid5);