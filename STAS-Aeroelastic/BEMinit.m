function x0 = BEMinit (psiFlag,                         ...
                       Viguess,Tar,Trg,Vg,wg,zr,ch,Lel, ...
                       aoast,rho,A,Dia,azi,omega)
%
% Gives a steady-state estimate of the states based on a Vi guess.
%
% Version:        Changes:
% --------        -------------
% 15.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 15.12.2017      
%
% Inputs:
% -------
% psiFlag         : Set to 1 to implement the dynamic wake in multi-blade
%                   coordinates.  Set to zero for blade-by-blade. 
%                   Regardless, dxdt is reported as blade-by-blade.
%                   If psiFlag = 1, then it is assumed that the elements
%                   are packed as Bl1 e1 e2 e3 ... eN, Bl2 e1 e2 ... for
%                   three blades.
% Viguess         : Virz and Virt for each element.
% Tar             : 3-by-3*Nel matrix of aero-to-rotorplane transforms.
% Trg             : Transform from rotorplane to global coordinates.
% Vg              : 3*Nel vector, incoming windspeed x,y,z in global
%                   coordinates.
% wg              : 3*Nel vector, structural motion, including rotor
%                   speed, x,y,z in global coordinates.
% zr              : 2*Nel vector, x,y projections of element positions
%                   onto the rotorplane.
% ch              : Airfoil chord length.
% rho             : Air density.
% A               : Area of the projected element on the rotorplane.
% Dia             : Diameter of the projected outer node on the rotorplane.
% omega           : Rotor speed.
%
% Outputs:
% --------
% x0              : Guess for initial x.

%{
'----------BEMinit-----------'
psiFlag
Viguess
Tar
Trg
Vg
wg
zr
ch
rho
A
Dia
omega
%}

Nel = size(Viguess,1)/2;
A0 = spalloc (7*Nel,7*Nel,12*Nel);
b0 = zeros (7*Nel,1);

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
Virxyz(i3a) = -Viguess(i2b).*spe;
Virxyz(i3b) =  Viguess(i2b).*cpe;
Virxyz(i3c) =  Viguess(i2a);

r   = sqrt(zr(i2a).^2 + zr(i2b).^2);

% Compute the local velocity.
Ua = zeros(3*Nel,1);
Ur = zeros(3*Nel,1);
Tgr = Trg.';
for iel = 1:Nel
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
% band.
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

for iel = 1:Nel
   j7 = 7*(iel-1);
   A0(j7+1,j7+1) = -tm1(iel);
   A0(j7+1,j7+2) =  tm1(iel)*K1(iel);
   A0(j7+1,j7+3) =  tm1(iel)*K2(iel);
   A0(j7+2,j7+3) =  1;
   A0(j7+3,j7+2) =  A32(iel);
   A0(j7+3,j7+3) =  A33(iel);
   b0(j7+1) = tm1(iel)*K3*aq(iel);
   b0(j7+3) = aq(iel);
end

% Wind velocity.
Vrxyz = zeros(3*Nel,1);
for iel = 1:Nel
   j3 = 3*(iel-1);
   Vrxyz(j3+[1:3]) = Tgr*(Vg(j3+[1:3]) - wg(j3+[1:3]));
end
Vrzts = zeros(3*Nel,1);
Vrzts(i3a) =  Vrxyz(i3c);
Vrzts(i3b) = -Vrxyz(i3a).*spe + Vrxyz(i3b).*cpe + r.*omega;  % Cancel r*omega in wg.
Vrzts(i3c) =  Vrxyz(i3a).*cpe + Vrxyz(i3b).*spe;

% Dynamic wake.
Vmag = sqrt(Vrzts(i3a).^2 + Vrzts(i3b).^2 + Vrzts(i3c).^2);
aa   = -Viguess(i2a)./Vrzts(i3a);

if (psiFlag == 0)

   tau1 = (1.1./(1 - 0.3*aa)).*(Dia./(2*Vmag));
   tau2 = (0.39 - 0.26*(2*r./Dia).^2).*tau1;
   it1  = 1./tau1;
   it2  = 1./tau2;

   for iel = 1:Nel

      j7 = 7*(iel-1);
      j2 = 2*(iel-1);

      A0(j7+4,j7+4) = -it1(iel);
      A0(j7+6,j7+4) =  it2(iel);
      A0(j7+6,j7+6) = -it2(iel);
      b0(j7+4)      =  0.4*it1(iel)*Viguess(j2+1);
      b0(j7+6)      =  0.6*it1(iel)*Viguess(j2+1);

      A0(j7+5,j7+5) = -it1(iel);
      A0(j7+7,j7+5) =  it2(iel);
      A0(j7+7,j7+7) = -it2(iel);
      b0(j7+5)      =  0.4*it1(iel)*Viguess(j2+2);
      b0(j7+7)      =  0.6*it1(iel)*Viguess(j2+2);

   end

elseif (psiFlag == 1)

   tau1 = (1.1./(1 - 0.3*aa)).*(Dia./(2*Vmag));
   tau2 = (0.39 - 0.26*(2*r./Dia).^2).*tau1;
   tau1 = mean(tau1)*ones(size(tau1,1),1);
   tau2 = mean(tau2)*ones(size(tau2,1),1);
   it1  = 1./tau1;
   it2  = 1./tau2;

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
      
      A0(i7d(ind),i7d(ind)) = psiB*(A1*Bpsi - omega*dBpsi);
      A0(i7e(ind),i7e(ind)) = psiB*(A1*Bpsi - omega*dBpsi);
      A0(i7f(ind),i7d(ind)) = psiB*A2*Bpsi;
      A0(i7f(ind),i7f(ind)) = psiB*(A3*Bpsi - omega*dBpsi);
      A0(i7g(ind),i7e(ind)) = psiB*A2*Bpsi;
      A0(i7g(ind),i7g(ind)) = psiB*(A3*Bpsi - omega*dBpsi);
      b0(i7d(ind))          = psiB*B1*Bpsi*Viguess(i2a(ind));
      b0(i7e(ind))          = psiB*B1*Bpsi*Viguess(i2b(ind));
      b0(i7f(ind))          = psiB*B2*Bpsi*Viguess(i2a(ind));
      b0(i7g(ind))          = psiB*B2*Bpsi*Viguess(i2b(ind));

   end

end

x0 = -A0\b0;

%{
psie
Ua
Ur
Uxy
aq
tau
Vrxyz
Vrzts
Vmag
aa
tau1
tau2
A0
b0
x0
'--------end BEMinit---------'
%}
