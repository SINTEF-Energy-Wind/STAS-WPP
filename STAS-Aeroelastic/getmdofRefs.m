function [imdofs,Nmd] = getmdofRefs (s)

imdofs(1) = 0;
imdofs(2) = 6 + s.foundation.Nmod;
imdofs(3) = imdofs(2) + s.tower.Nmod;
imdofs(4) = imdofs(3) + s.nacelle.Nmod;
imdofs(5) = 0;
imdofs(6) = imdofs(4) + s.driveshaft.Nmod;
imdofs(7) = imdofs(6) + s.blade(1).Nmod;
imdofs(8) = imdofs(7) + s.blade(2).Nmod;

Nmd = 12 + s.foundation.Nmod + s.tower.Nmod + s.nacelle.Nmod ...
    + s.driveshaft.Nmod + s.blade(1).Nmod + s.blade(2).Nmod  ...
    + s.blade(3).Nmod;