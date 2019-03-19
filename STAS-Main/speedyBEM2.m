function [Vir,Vmag,f,aoa,phi,Cl,Cd,dCl,dCd,Fld,Fas,Fap] =           ...
                            speedyBEM2 (Neb,Vinf,omega,psi,beta,xi, ...
                                        Dia,r,Lel,chord,xpc,        ...
                                        aoas,kfoils,foilwt,         ...
                                        density,viscosity)
%
% This subroutine performs a very simple static BEM calculation
% that does not include yawed flow.
%
% Version:        Changes:
% --------        -------------
% 10.06.2016      Rewritten based on speedyBEM.m, employing fsolve instead
%                 of my relaxation method.
% 30.06.2017      Adapted for complex step.
% 08.08.2017      Updated airfoil coefficient computation to employ a new,
%                 easier spline format.  Got rid of global variables.
%
% Version:        Verification:
% --------        -------------
% 10.06.2016      
% 30.06.2017      
% 08.08.2017      
%
% Inputs:
% -------
% Vinf            : Incoming axial velocity at each element.
% omega           : Rotational speed.
% psi             : Driveshaft azimuth angle.
% beta            : Vector length 3, pitch angle for each blade.
% xi              : Twist angles (rad).
% Dia             : Rotor diameter.
% r               : Element radii from axis of rotation.
% Lel             : Element lengths.
% chord           : Element chords.
% xpc             : Distance from LE to pitch axis, fraction of chord.
% aoas            : Angles-of-attack in data file, radians.
% kfoils          : Spline coefficients describing airfoil coefficients.
%                   Generate these with runACS.m.
% foilwt          : An Nfoils-by-NN array containing the weight to
%                   apply to each airfoil type's coefficients, for
%                   each blade element.
% density, viscosity
%
% Outputs:
% --------
% Vir             : Axial and tangential induced velocity in
%                   rotorplane coordinates.
% Vmag            : Local velocity at each element.
% f               : Prandtl factor.
% aoa             : Angle-of-attack (rad).
% phi             : Inflow angle (rad).
% Cl,Cd           : Lift and drag coefficients.
% dCl,dCd         : dCl/daoa, dCd/daoa.
% Fld             : Vector of (Fl,Fd) lift/drag force pairs.
% Fas             : Airfoil forces and moments in section
%                   coordinates.
% Fap             : Airfoil forces and moments in pitch
%                   coordinates.

NN = size(Vinf)(1);

Vir = zeros(2*NN,1);
Fld = zeros(2*NN,1);
Fas = zeros(6*NN,1);
Fap = zeros(6*NN,1);

xib = xi + [beta(1)*ones(Neb,1);beta(2)*ones(Neb,1);beta(3)*ones(Neb,1)];

if (real(omega) == 0)

   % Compute rotor aero loads, no induced velocity, so BEM is not needed.
   phi = pi/2*ones(NN,1);
   aoa = phi - xib;
   V = absc(Vinf);
   C = zeros(6,NN);
   for iel = 1:NN
      C(:,iel) = airfoilCoefficientsSpline (aoas,kfoils,foilwt(:,iel),aoa(iel));
   end
   Cl = C(1,:).';
   Cd = C(2,:).';

   Viz = zeros(NN,1);
   Vit = zeros(NN,1);
   f = ones(NN,1);

else

   Viz = zeros(NN,1);
   Vit = zeros(NN,1);

   for iel = 1:NN

      % Initial guess for induced velocities.
      Viz0 = -Vinf(iel)/3;
      Vit0 = 0;
      vinit = [Viz0;Vit0];

      [v,fvec,info,output,fjac] = ...
         fsolve (@(Vi) fbem (Vi,Vinf(iel),Dia,r(iel),omega,chord(iel), ...
                             aoas,kfoils,foilwt(:,iel),xib(iel)), vinit);

      Viz(iel) = v(1);
      Vit(iel) = v(2);

   end

end

Vz = Vinf + Viz;
Vt = r(1:NN).*omega + Vit;
Vmag = sqrt(Vz.^2 + Vt.^2);

% Compute final forces and induced velocity.
val = 0.5d0*density*chord(1:NN).*Lel(1:NN).*Vmag.^2;
phi = atan2c(Vz,Vt);
cp = cos(phi);
sp = sin(phi);
aoa = phi - xib(1:NN);
C = zeros(6,NN);
for iel = 1:NN
   C(:,iel) = airfoilCoefficientsSpline (aoas,kfoils,foilwt(:,iel),aoa(iel));
end
sa = sin(aoa);
ca = cos(aoa);
sx = sin(xi(1:NN));
cx = cos(xi(1:NN));
Cl = C(1,:).';
Cd = C(2,:).';
Cm = C(3,:).';
dCl = C(4,:).';
dCd = C(5,:).';

f = zeros(NN,1);
ind = (real(sp) < 1e-3 & real(sp) >= 0) | (real(sp) < 0);
nind = ~ind;
f(ind) = 1;
f(nind) = maxc((2/pi)*acos(exp(-3*(0.5*Dia - r(nind)) ...
        ./(2*r(nind).*sp(nind)))),0.1*ones(size(r(nind),1),1));

%ifx = [1:6:6*NN-5].';
ify = [2:6:6*NN-4].';
ifz = [3:6:6*NN-3].';
imx = [4:6:6*NN-2].';
%imy = [5:6:6*NN-1].';
%imz = [6:6:6*NN].';

Fas(ify) = val.*(Cl.*sa - Cd.*ca);
Fas(ifz) = val.*(Cl.*ca + Cd.*sa);
Fas(imx) = val.*chord(1:NN).*Cm ...
         + Fas(ifz).*(xpc(1:NN) - 0.25).*chord(1:NN);

Fap(ify) =  Fas(ify).*cx + Fas(ifz).*sx;
Fap(ifz) = -Fas(ify).*sx + Fas(ifz).*cx;
Fap(imx) =  Fas(imx);

iz = [1:2:2*NN-1].';
it = [2:2:2*NN].';

Fld(iz) = val.*Cl;
Fld(it) = val.*Cd;

Vir(iz) = Viz;
Vir(it) = Vit;




