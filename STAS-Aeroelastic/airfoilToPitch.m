function [Fp,Dy] = airfoilToPitch (Fa,TsB,dTsB,Tas)
%
% The transform of a vector from airfoil to pitch coordinates.
%
%   States:           y vector:         u vector:
%                     qn1,2     1:12
%                     Fa       13:18
%
% Version:        Changes:
% --------        -------------
% 24.02.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 24.02.2018      Double-checked code.
%
% Inputs:
% -------
% Fa              : Forces in airfoil coordinates.
% TsB,dTsB        : Section-to-body transform and its derivatives.
% Tas             : Airfoil-to-section transform, a constant matrix
%                   for a given element.
%
% Outputs:
% --------
% Fp

Dy = spalloc(6,18);

Fp = zeros(6,1);

TT = TsB*Tas;
T6 = [TT zeros(3,3);zeros(3,3) TT];
Fp = T6*Fa;

Dy(:,13:18) = T6;

for jj = 1:12

   jc3 = 3*(jj-1);

   dTT = dTsB(:,jc3+[1:3])*Tas;
   Dy(:,jj) = [dTT zeros(3,3);zeros(3,3) dTT]*Fa;

end
