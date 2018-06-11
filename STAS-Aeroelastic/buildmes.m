function mes = buildmes (rhos,L)
%
% Ordering of DOFs is
% wxs1,wys1,wzs1,txs1,tys1,tzs1,wxs2,wys2,wzs2,txs2,tys2,tzs2
%
% rhos is a 6-by-6 matrix, rhos = int{rho (R^T)R}dA}.  See
% Merz KO.  Enhanced finite beam elements for wind turbine blades.
% Memo AN 16.12.18, SINTEF Energy Research, 2016.
%
% Version:        Changes:
% --------        -------------
% 30.10.2014      Original code.
% 15.04.2016      Modified section properties to matrix format.
% 30.06.2017      Adapted for complex step.
%
% Version:        Verification:
% --------        -------------
% 30.10.2014      Verified natural frequencies in bending, axial,
%                 and torsion against analytical solutions, using
%                 one and two elements.
% 15.04.2016      Reproduces results of previous version.
% 30.06.2017      
%

% Define coefficients on S polynomial.

S3 = spalloc(6,12,8);
S2 = spalloc(6,12,16);
S1 = spalloc(6,12,14);
S0 = spalloc(6,12,6);

S3(2,2)  = 2/(L^3);
S3(2,6)  = 1/(L^2);
S3(2,8)  = -2/(L^3);
S3(2,12) = 1/(L^2);
S3(3,3)  = S3(2,2);
S3(3,5)  = -S3(2,6);
S3(3,9)  = S3(2,8);
S3(3,11) = -S3(2,12);

S2(2,2)  = -3/(L^2);
S2(2,6)  = -2/L;
S2(2,8)  = 3/(L^2);
S2(2,12) = -1/L;
S2(3,3)  = S2(2,2);
S2(3,5)  = -S2(2,6);
S2(3,9)  = S2(2,8);
S2(3,11) = -S2(2,12);
S2(5,3)  = -6/(L^3);
S2(5,5)  = 3/(L^2);
S2(5,9)  = -S2(5,3);
S2(5,11) = 3/(L^2);
S2(6,2)  = -S2(5,3);
S2(6,6)  = S2(5,5);
S2(6,8)  = S2(5,3);
S2(6,12) = S2(5,11);

S1(1,1)  = -1/L;
S1(1,7)  = 1/L;
S1(2,6)  = 1;
S1(3,5)  = -S1(2,6);
S1(4,4)  = S1(1,1);
S1(4,10) = S1(1,7);
S1(5,3)  = 6/(L^2);
S1(5,5)  = -4/L;
S1(5,9)  = -S1(5,3);
S1(5,11) = -2/L;
S1(6,2)  = -S1(5,3);
S1(6,6)  = S1(5,5);
S1(6,8)  = S1(5,3);
S1(6,12) = S1(5,11);

S0(1,1)  = 1;
S0(2,2)  = 1;
S0(3,3)  = S0(2,2);
S0(4,4)  = S0(1,1);
S0(5,5)  = 1;
S0(6,6)  = S0(5,5);

% Multiply to get the new coefficients.
SRS6 = (S3.')*rhos*S3;
SRS5 = (S2.')*rhos*S3 + (S3.')*rhos*S2;
SRS4 = (S2.')*rhos*S2 + (S3.')*rhos*S1 + (S1.')*rhos*S3;
SRS3 = (S2.')*rhos*S1 + (S1.')*rhos*S2 + (S3.')*rhos*S0 + (S0.')*rhos*S3;
SRS2 = (S2.')*rhos*S0 + (S0.')*rhos*S2 + (S1.')*rhos*S1;
SRS1 = (S1.')*rhos*S0 + (S0.')*rhos*S1;
SRS0 = (S0.')*rhos*S0;

% Integrate and evaluate at L.
mes = L*(SRS0 + L*((SRS1/2) + L*((SRS2/3) + L*((SRS3/4) ...
    + L*((SRS4/5) + L*((SRS5/6) + L*(SRS6/7)))))));