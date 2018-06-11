function Pin = assemblePin (s)
% Version:        Changes:
% --------        -------------
% 01.02.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 01.02.2018      undeformedPosition uses the output Pin.

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

Pin = zeros(Ndof,1);

idof = 0;

Pin(idof+[1:s.foundation.Ndof]) = s.foundation.Pn_B;
idof = idof + s.foundation.Ndof;

Pin(idof+[1:s.tower.Ndof]) = s.tower.Pn_B;
idof = idof + s.tower.Ndof;

Pin(idof+[1:s.nacelle.Ndof]) = s.nacelle.Pn_B;
idof = idof + s.nacelle.Ndof;

Pin(idof+[1:s.driveshaft.Ndof]) = s.driveshaft.Pn_B;
idof = idof + s.driveshaft.Ndof;

Pin(idof+[1:s.blade(1).Ndof]) = s.blade(1).Pn_B;
idof = idof + s.blade(1).Ndof;

Pin(idof+[1:s.blade(2).Ndof]) = s.blade(2).Pn_B;
idof = idof + s.blade(2).Ndof;

Pin(idof+[1:s.blade(3).Ndof]) = s.blade(3).Pn_B;
idof = idof + s.blade(3).Ndof;
