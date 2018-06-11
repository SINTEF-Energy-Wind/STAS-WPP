function Cl = getCL (aoas,kfoils,foilwt,aoa)
% Used with fzero to find aoaz.

C = airfoilCoefficientsSpline (aoas,kfoils,foilwt,aoa);
Cl = C(1);
