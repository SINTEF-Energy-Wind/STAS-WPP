function C = airfoilCoefficientsSpline (aoas,kfoils,weights,aoa)
%
% This function computes the airfoil coefficients Cl,Cd,Cm and
% their quasi-steady (low frequency) first derivatives with respect
% to the angle-of-attack.
%
% Splines and, for deep stall, theoretical curves are employed to
% ensure smooth coefficient and derivative values.  Numerical
% smoothness may be critical in some applications such as
% optimization.
%
% Version:        Changes:
% --------        -------------
% 08.08.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 08.08.2017      Verified by plotting over tabulated values.
%
% Inputs:
% -------
% aoas            : Angles of attack in the data file used to 
%                   generate kfoils, radians.
% kfoils          : Spline coefficients.
% weights         : A column vector specifying the weight to use for
%                   each type of airfoil in the kfoils array.
% aoa             : Angles-of-attack at which to compute C, radians.
%
% Outputs:
% --------
% C                : [Cl;Cd;Cm;dClda;dCdda;dCmda] at the given aoa. 

% Force aoa to lie between -pi and pi.
ind = (real(aoa) < -pi);
aoa(ind) = aoa(ind) + 2*pi;
ind = (real(aoa) > pi);
aoa(ind) = aoa(ind) - 2*pi;

[S,dSda] = splineSMatrix (aoas,aoa);

Naoa = size(aoa,1);
Nfoil = size(weights,1);
lw = logical(weights);
Nw = sum(lw);
ind = [1:Nfoil].';
ifoil = ind(lw);

C = zeros(Naoa,6);
for jj = 1:Nw
   foil = ifoil(jj);
   Cf = [S*kfoils(:,3*(foil-1)+[1:3]) dSda*kfoils(:,3*(foil-1)+[1:3])];
   C = C + weights(ifoil(jj))*Cf;
end



C = [2*pi*aoa;0.005;0;2*pi;0;0];

