function [dxdt,y] = BEMNL (psiFlag,                       ...
                           x,t,Tar,Trg,Vg,wg,zr,ch,Lel,   ...
                           aoas,kfoils,foilwt,aoaz,aoast, ...
                           xas,yas,rho,A,Dia,azi,omega)
%
% State-space equations for the blade element momentum method.  Vectorized
% to the extent possible for multiple element inputs.
%
%   States:           y vector:
%   ad        1       Fa        1:6
%   a1,a2    2:3      
%   Vih z,t  4:5      
%   Vi z,t   6:7
%
% Version:        Changes:
% --------        -------------
% 12.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 12.12.2017      
%
% Inputs:
% -------
% psiFlag         : Set to 1 to implement the dynamic wake in multi-blade
%                   coordinates.  Set to zero for blade-by-blade. 
%                   Regardless, dxdt is reported as blade-by-blade.
%                   If psiFlag = 1, then it is assumed that the elements
%                   are packed as Bl1 e1 e2 e3 ... eN, Bl2 e1 e2 ... for
%                   three blades.
% x               : The vector of states.  For each element, these are
%                   ad, a1, a2, Virhz, Virht, Virz, and Virt.
% Tar             : 3-by-3*Nel matrix of aero-to-rotorplane transforms.
% Trg             : Transform from rotorplane to global coordinates.
% Vg              : 3*Nel vector, incoming windspeed x,y,z in global
%                   coordinates.
% wg              : 3*Nel vector, structural motion, including rotor
%                   speed, x,y,z in global coordinates.
% zr              : 2*Nel vector, x,y projections of element positions
%                   onto the rotorplane.
% ch              : Airfoil chord length.
% aoas,kfoils     : Splined airfoil coefficient tables.
% foilwt          : Nfoil-by-Nel table.  Weights to use when computing
%                   airfoil coefficients from splined tables.
% aoaz            : Zero-lift angles-of-attack for each element. Should
%                   be computed precisely from the airfoil tables.
% aoast           : 2*Nel vector, containing deep-stall angles-of-attack
%                   for each element.  Alternating positive, negative.
% xas,yas         : X^a and Y^a coordinates of the reference section
%                   coordinate system.
% rho             : Air density.
% A               : Area of the projected element on the rotorplane.
% Dia             : Diameter of the projected outer node on the rotorplane.
% azi             : A measure of the rotor azimuth, such that 
%                   dpsi/dt = omega.
% omega           : Rotor speed.
%
% Outputs:
% --------
% dxdt            : Rate of change of states x.
% y               : Vector of outputs which are to be passed to other
%                   modules.

%tic
%del = sqrt(eps);

dxdt = zeros(size(x,1),1);

Nel = size(x,1)/7;

i2a = [1:2:2*Nel-1].';
i2b = [2:2:2*Nel].';
i3a = [1:3:3*Nel-2].';
i3b = [2:3:3*Nel-1].';
i3c = [3:3:3*Nel].';
i6a = [1:6:6*Nel-5].';
i6b = [2:6:6*Nel-4].';
i6c = [3:6:6*Nel-3].';
i6d = [4:6:6*Nel-2].';
i6e = [5:6:6*Nel-1].';
i6f = [6:6:6*Nel].';
i7a = [1:7:7*Nel-6].';
i7b = [2:7:7*Nel-5].';
i7c = [3:7:7*Nel-4].';
i7d = [4:7:7*Nel-3].';
i7e = [5:7:7*Nel-2].';
i7f = [6:7:7*Nel-1].';
i7g = [7:7:7*Nel].';

psie = atan2c (zr(i2b),zr(i2a));
cpe = cos(psie);
spe = sin(psie);
Virxyz = zeros(3*Nel,1);
Virxyz(i3a) = -x(i7g).*spe;
Virxyz(i3b) =  x(i7g).*cpe;
Virxyz(i3c) =  x(i7f);

% Compute the local velocity.
Ua = zeros(3*Nel,1);
Ur = zeros(3*Nel,1);
Tgr = Trg.';
for iel = 1:Nel  % A necessary evil, due to the transform matrix.
   j3 = 3*(iel-1);
   Tra = Tar(:,j3+[1:3]).';
   Ur(j3+[1:3]) = Virxyz(j3+[1:3]) + Tgr*(Vg(j3+[1:3]) - wg(j3+[1:3]));
   Ua(j3+[1:3]) = Tra*Ur(j3+[1:3]);
end
Uxy = sqrt(Ua(i3a).^2 + Ua(i3b).^2);  % Neglect spanwise component for airfoils.

% Compute the qs angle-of-attack.
aq = atan2c (Ua(i3b),Ua(i3a));
caq = cos(aq);
saq = sin(aq);

% Dynamic aoa parameters.
uc  = 2*Uxy./ch;
A1  = 0.165;
A2  = 0.335;
b1  = 0.0455;
b2  = 0.3;
K1  = (A1 + A2)*b1*b2*uc.^2;
K2  = (A1*b1 + A2*b2)*uc;
K3  = 1 - A1 - A2;
A32 = -b1*b2*uc.^2;
A33 = -(b1+b2)*uc;
tau = 8.6./uc;

%{
% Special logic: tau is decreased at high angles-of-attack such that
% the dynamic and quasi-steady aoa's match over the relevant frequency
% band.  [This causes problems for the matching linearization, and it's
% probably not worth the trouble for the sorts of cases I'll be looking
% at.]
aoastp = aoast(i2a);
aoastn = aoast(i2b);
ind = (real(aq) > aoastp);
sind = sum(ind);
tau(ind) = minc(maxc((1 - ((aq(ind) - aoastp(ind))/0.1)*0.8), ...
                     0.2*ones(sind,1)), ones(sind,1)).*tau(ind);
ind = (real(aq) < aoastn);
sind = sum(ind);
tau(ind) = minc(maxc((1 + ((aq(ind) - aoastn(ind))/0.1)*0.8), ...
                     0.2*ones(sind,1)), ones(sind,1)).*tau(ind);
%}

tm1 = 1./tau;
dxdt(i7a) = -tm1.*x(i7a) + tm1.*K1.*x(i7b) + tm1.*K2.*x(i7c) + tm1.*K3.*aq;
dxdt(i7b) = x(i7c);
dxdt(i7c) = A32.*x(i7b) + A33.*x(i7c) + aq;

% Baseline airfoil coefficients: Cl,Cd,Cm,dCl/da,dCd/da,dCm/da.
C  = zeros(Nel,6);
Cq = zeros(Nel,6);
Cl = zeros(Nel,1);
for iel = 1:Nel
   i7 = 7*(iel-1);
   C(iel,:)  = airfoilCoefficientsSpline (aoas,kfoils,foilwt(:,iel),x(i7+1));
   Cq(iel,:) = airfoilCoefficientsSpline (aoas,kfoils,foilwt(:,iel),aq(iel));
   if (absc(x(i7+1) - aoaz(iel)) < eps^(0.25))
      Cl(iel) = C(iel,4)*(aq(iel) - x(i7+1));
   else
      Cl(iel) = (1 + (aq(iel) - x(i7+1))/(x(i7+1) - aoaz(iel)))*C(iel,1);
   end
end

Cd = Cq(:,2);
Cm = Cq(:,3);

val = 0.5*rho*ch.*Lel.*(Uxy.^2);
Fa = zeros(6*Nel,1);
Fa(i6a) = val.*(-Cl.*saq + Cd.*caq);
Fa(i6b) = val.*(Cl.*caq + Cd.*saq);
Fa(i6f) = -val.*ch.*Cm + Fa(i6a).*yas - Fa(i6b).*xas;

Frxyz = zeros(3*Nel,1);
for iel = 1:Nel
   j3 = 3*(iel-1);
   j6 = 6*(iel-1);
   Frxyz(j3+[1:3]) = Tar(:,j3+[1:3])*Fa(j6+[1:3]);
end

Frzt = zeros(2*Nel,1);
Frzt(i2a) = Frxyz(i3c);
Frzt(i2b) = -Frxyz(i3a).*spe + Frxyz(i3b).*cpe;

% Momentum balance.  First the Prandtl factor.
Uz  =  Ur(i3c);
Ut  = -Ur(i3a).*spe + Ur(i3b).*cpe;
r   = sqrt(zr(i2a).^2 + zr(i2b).^2);

pr  = prandtl (Uz,Ut,r,Dia);

% Quasi-steady induced velocity.
Vrxyz = zeros(3*Nel,1);
for iel = 1:Nel
   j3 = 3*(iel-1);
   Vrxyz(j3+[1:3]) = Tgr*(Vg(j3+[1:3]) - wg(j3+[1:3]));
end
Vrzts = zeros(3*Nel,1);
Vrzts(i3a) =  Vrxyz(i3c);
Vrzts(i3b) = -Vrxyz(i3a).*spe + Vrxyz(i3b).*cpe + r.*omega;  % Cancel r*omega.
Vrzts(i3c) =  Vrxyz(i3a).*cpe + Vrxyz(i3b).*spe;
Vir = zeros(2*Nel,1);
Vir(i2a) = x(i7f);
Vir(i2b) = x(i7g);

Viq = ViqsNL (Frzt,Vrzts,Vir,r,Dia,rho,A,pr);

% Dynamic wake.
Vmag = sqrt(Vrzts(i3a).^2 + Vrzts(i3b).^2 + Vrzts(i3c).^2);
aa   = -Vir(i2a)./Vrzts(i3a);

if (psiFlag == 0)

   tau1 = (1.1./(1 - 0.3*aa)).*(Dia./(2*Vmag));
   tau2 = (0.39 - 0.26*(2*r./Dia).^2).*tau1;
   it1  = 1./tau1;
   it2  = 1./tau2;

   dxdt(i7d) = -it1.*x(i7d) + 0.4*it1.*Viq(i2a);
   dxdt(i7e) = -it1.*x(i7e) + 0.4*it1.*Viq(i2b);
   dxdt(i7f) =  it2.*x(i7d) - it2.*x(i7f) + 0.6*it2.*Viq(i2a);
   dxdt(i7g) =  it2.*x(i7e) - it2.*x(i7g) + 0.6*it2.*Viq(i2b);

elseif (psiFlag == 1)

   tau1 = (1.1./(1 - 0.3*aa)).*(Dia./(2*Vmag));
   tau2 = (0.39 - 0.26*(2*r./Dia).^2).*tau1;
   tm1 = mean(tau1)*ones(size(tau1,1),1);
   tm2 = mean(tau2)*ones(size(tau2,1),1);
   it1  = 1./tm1;
   it2  = 1./tm2;

   caz = [cos(azi) cos(azi + 2*pi/3) cos(azi + 4*pi/3)].';
   saz = [sin(azi) sin(azi + 2*pi/3) sin(azi + 4*pi/3)].';

   psiB  = [1 caz(1) saz(1); ...
            1 caz(2) saz(2); ...
            1 caz(3) saz(3)];
   Bpsi  = (1/3)*[   1         1         1    ; ...
                  2*caz(1)  2*caz(2)  2*caz(3); ...
                  2*saz(1)  2*saz(2)  2*saz(3)];
   dBpsi = (1/3)*[   0         0         0    ; ...
                 -2*saz(1) -2*saz(2) -2*saz(3); ...
                  2*caz(1)  2*caz(2)  2*caz(3)];

   Nb = 3;
   Neb = Nel/Nb;
   for iel = 1:Neb
      ind = [iel Neb+iel 2*Neb+iel].';

      A1 = diag(-it1(ind));
      A2 = diag(it2(ind));
      A3 = diag(-it2(ind));
      B1 = 0.4*diag(it1(ind));
      B2 = 0.6*diag(it2(ind));
      dxdt(i7d(ind)) = psiB*(A1*Bpsi - omega*dBpsi)*x(i7d(ind)) ...
                     + psiB*B1*Bpsi*Viq(i2a);
      dxdt(i7e(ind)) = psiB*(A1*Bpsi - omega*dBpsi)*x(i7e(ind)) ...
                     + psiB*B1*Bpsi*Viq(i2b);
      dxdt(i7f(ind)) = psiB*A2*Bpsi*x(i7d(ind))                 ...
                     + psiB*(A3*Bpsi - omega*dBpsi)*x(i7f(ind)) ...
                     + psiB*B2*Bpsi*Viq(i2a);
      dxdt(i7g(ind)) = psiB*A2*Bpsi*x(i7e(ind))                 ...
                     + psiB*(A3*Bpsi - omega*dBpsi)*x(i7g(ind)) ...
                     + psiB*B2*Bpsi*Viq(i2b);

   end

end

y = Fa;
