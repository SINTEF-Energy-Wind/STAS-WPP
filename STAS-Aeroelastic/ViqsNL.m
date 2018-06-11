function Viq = ViqsNL (Fr,Vr,Vir,r,D,rho,A,f)
%
% Nonlinear equation for the quasi-steady induced velocity, including
% Glauert's formula.  Vectorized if there are multiple blade elements.
%
% Version:        Changes:
% --------        -------------
% 11.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 11.12.2017      
%
% Inputs:
% -------
% Fr              : Blade element forces F z,t, rotor coordinates.
% Vr              : Incoming windspeed V z,t,s, rotor coordinates.
% Vir             : Induced velocity Vi z,t, rotor coordinates.
% r,D             : Radius, rotor outer diameter (projected onto the
%                   rotor plane, if the blades are deformed).
% rho             : Air density.
% A               : Swept area associated with the blade element.
% f               : Tip loss function.
% 
% Outputs:
% --------
% Viq             : Viq z,t.

Nel = size(Fr,1)/2;

Viqa = zeros(2*Nel,1);  % An "axial" value, before the yaw correction.
Viq  = zeros(2*Nel,1);

od = [1:2:2*Nel-1].';
ev = [2:2:2*Nel].';
jz = [1:3:3*Nel-2].';
jt = [2:3:3*Nel-1].';
js = [3:3:3*Nel].';

W = sqrt((Vr(jz) + f.*Vir(od)).^2 + Vr(jt).^2 + Vr(js).^2);
Vrmag = sqrt(Vr(jz).^2 + Vr(jt).^2 + Vr(js).^2);

Viqa(ev) = -Fr(ev)./(2.*rho.*A.*f.*W);
Viqa(od) = -Fr(od)./(2.*rho.*A.*f.*W);
a2 = 1;
CT2 = 1.82;
a1 = 1 - 0.5*sqrt(CT2);
CT1 = 4*a1*(1 - a1);
ind = (real(Viqa(od)) < real(-(a1./f).*Vr(jz)));

if (sum(ind) > 0)
   Viqaz = Viqa(od);
   Vrz = Vr(jz);
   Frz = Fr(od);
   Viqaz(ind) = -(Vrz(ind)./f(ind))                                         ...
             .* ((Frz(ind)./(0.5*(CT2 - CT1).*rho.*A(ind).*(Vrmag(ind).^2)) ...
              -   CT1/(CT2 - CT1))*(a2 - a1) + a1);
   Viqa(od) = Viqaz;
end

mag = sqrt(Vr(jt).^2 + Vr(js).^2);
Vrs = Vr(js);
ind = real(mag) > sqrt(eps);
cp0 = zeros(Nel,1);
cp0(ind) = Vrs(ind)./mag(ind);
tanv = tan(0.5*acos((Vr(jz) + f.*Vir(od))./W));
fyaw = 1 + (2*r./D).*tanv.*cp0;

Viq(od) = fyaw.*Viqa(od);
Viq(ev) = fyaw.*Viqa(ev);

%-Fr(od)./(2*rho*A.*f.*sqrt((Vr(jz) + f.*Viq(od)).^2 + Vr(jt).^2 + Vr(js).^2)) - Viq(od)
%-Fr(ev)./(2*rho*A.*f.*sqrt((Vr(jz) + f.*Viq(od)).^2 + Vr(jt).^2 + Vr(js).^2)) - Viq(ev)

