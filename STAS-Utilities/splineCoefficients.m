function [A,B] = splineCoefficients (r)
%
% This function formulates the matrices in the spline fitting problem
%   ^        ^  ^
% A(r)*k = B(r)*z,                               ^ ^
% where k are the cubic spline coefficients and (r,z) are the control
% points.
%
% The function is intentionally formulated inefficiently, from a 
% computational standpoint, such that it is easy to formulate in the
% A*k = B*z format.  The formulation is based on
% z = k0 + k1*s + k2*s^2 + k3*s^3,      0<=s<=1,
% s := (x-x_{i-1})/h_i,   h_i := x_i - x_{i-1}.
%
% This means that k0 is the value of z at x = x_{i-1}, and
% k0 + k1 + k2 + k3 is the value of z at x = x_i.
%
% Version:        Changes:
% --------        -------------
% 03.12.2016      Original code.
% 04.07.2017      Adapted for complex step derivatives.
%
% Version:        Verification:
% --------        -------------
% 03.12.2016      Memo AN 16.12.70.
% 04.07.2017      
%
% Inputs:
% -------
% r               : A list of x coordinates of the control points.
%
% Outputs:
% --------
% A,B             : Matrices satisfying the above equation.

Np = size(r,1);
Nk = 4*(Np-1);

h = r(2:Np) - r(1:Np-1);

ii = [1 1                                         ... % z''_0 = 0.5*z''_1
      [2:4:Nk-2]                                  ... % Match z_{i-1}.
      [3:4:Nk-1] [3:4:Nk-1] [3:4:Nk-1] [3:4:Nk-1] ... % Match z_i.
      [4:4:Nk-4] [4:4:Nk-4] [4:4:Nk-4] [4:4:Nk-4] ... % z'_i = z'_i.
      [5:4:Nk-3] [5:4:Nk-3] [5:4:Nk-3]            ... % z''_i = z''_i.
      Nk Nk].';                                       % z''_Np = 0.5*z''_{Np-1}

jj = [3 4                                         ...
      [1:4:Nk-3]                                  ...
      [1:4:Nk-3] [2:4:Nk-2] [3:4:Nk-1] [4:4:Nk]   ...
      [2:4:Nk-6] [3:4:Nk-5] [4:4:Nk-4] [6:4:Nk-2] ...
      [3:4:Nk-5] [4:4:Nk-4] [7:4:Nk-1]            ...
      Nk-1 Nk].';

ss = [1/(h(1)^2) -3/(h(1)^2)                                        ...
      ones(1,Np-1)                                                  ...
      ones(1,Np-1) ones(1,Np-1) ones(1,Np-1) ones(1,Np-1)           ...
      1./h(1:Np-2).' 2./h(1:Np-2).' 3./h(1:Np-2).' -1./h(2:Np-1).'  ...
      2./(h(1:Np-2).').^2 6./(h(1:Np-2).').^2 -2./(h(2:Np-1).').^2  ...
      1/h(Np-1)^2 6/h(Np-1)^2].';

A = sparse(ii,jj,ss,Nk,Nk);

ii = [[2:4:Nk-2] [3:4:Nk-1]].';                        % Match z_{i-1} and z_i.

jj = [[1:Np-1] [2:Np]].';

ss = [ones(1,Np-1) ones(1,Np-1)].';

B = sparse(ii,jj,ss,Nk,Np);

        