function f = fbem (Vi,VV,DD,ri,WW,chord,aoas,kfoils,foilwt,xib)
%
% Vi = [Viz Vit].'
%
% Version:        Changes:
% --------        -------------
% 30.06.2017      Adapted for complex step.
% 08.08.2017      Updated airfoil coefficient computation to employ a new,
%                 easier spline format.  Got rid of global variables.
%
% Version:        Verification:
% --------        -------------
% 30.06.2017      
% 08.08.2017      

f = zeros(2,1);

a2 = 1;
Ct2 = 1.82;
a1 = 1 - 0.5*sqrt(Ct2);
Ct1 = 4*a1*(1 - a1);

Vz = VV + Vi(1);
rw = ri*WW;
Vt = rw + Vi(2);
sol = 3*chord/(2*pi*ri);

phi = atan2c(Vz,Vt);
cp = cos(phi);
sp = sin(phi);
aoa = phi - xib;
Vmag = sqrt(Vz^2 + Vt^2);

C = airfoilCoefficientsSpline (aoas,kfoils,foilwt,aoa);

Cl = C(1);
Cd = C(2);

if ((real(sp) < 1e-3) && (real(sp) >= 0)) || (real(sp) < 0)
   Pr = 1;
else
   Pr = maxc((2/pi)*acos(exp(-3*(0.5*DD - ri) ...
      /      (2*ri*sp))),0.1);
end

val = 4*(VV + Pr*Vi(1))*Pr;

Ct = sol*((Vmag/VV)^2)*(Cl*cp + Cd*sp);

if (Ct > Ct1)
   f(1) = -Vi(1) - (VV/Pr)*(((Ct - Ct1)/(Ct2 - Ct1)) ...
        * (a2 - a1) + a1);
else
   f(1) = -Vi(1) + minc(-sol*(Vmag^2)*(Cl*cp + Cd*sp)/val,VV);
end

f(2) = -Vi(2) + minc(-sol*(Vmag^2)*(-Cl*sp + Cd*cp)/val,VV);

