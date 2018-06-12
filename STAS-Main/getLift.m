function Cl = getLift (aoas,kfoils,weights,aoa);

C     = airfoilCoefficientsSpline (aoas,kfoils,weights,aoa);
Cl    = C(:,1);


