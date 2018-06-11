function ke = buildkes (EEs,L)
%
% Ordering of DOFs is
% wxs1,wys1,wzs1,txs1,tys1,tzs1,wxs2,wys2,wzs2,txs2,tys2,tzs2
%
% EEs is a 6-by-6 matrix, EEs = int{(D^T)ED}dA}.  See
% Merz KO.  Enhanced finite beam elements for wind turbine blades.
% Memo AN 16.12.18, SINTEF Energy Research, 2016.
%
% Version:        Changes:
% --------        -------------
% 27.10.2014      Original code.
% 15.04.2016      Modified section properties to matrix format.
% 30.06.2017      Adapted for complex step.
%
% Version:        Verification:
% --------        -------------
% 27.10.2014      Verified natural frequencies in bending, axial,
%                 and torsion against analytical solutions, using
%                 one and two elements.
% 15.04.2016      Reproduces results of previous version.
% 30.06.2017      

% Define coefficients on S polynomial.

[S0,S1,S2,S3] = kscoeff (L);

% Compute the derivative.
dS2 = 3*S3;
dS1 = 2*S2;
dS0 = S1;

% Multiply to get the new coefficients:
SES4 = (dS2.')*EEs*dS2;
SES3 = (dS2.')*EEs*dS1 + (dS1.')*EEs*dS2;
SES2 = (dS2.')*EEs*dS0 + (dS0.')*EEs*dS2 ...
     + (dS1.')*EEs*dS1;
SES1 = (dS1.')*EEs*dS0 + (dS0.')*EEs*dS1;
SES0 = (dS0.')*EEs*dS0;

% Integrate and evaluate at L.
ke = L*(SES0 + L*((SES1/2) + L*((SES2/3) + L*((SES3/4) + L*(SES4/5)))));
