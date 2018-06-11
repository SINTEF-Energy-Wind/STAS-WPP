function [S,dSdx] = splineSMatrix (rhat,x)
%
% With splineCoefficients.m I compute a set of polynomial coefficients
% in vector form, such that (if I knew the zhat's) k = A\B*zhat.
% A point z can be related to the coefficients by
% z = [... 1 s s^2 s^3 ...] k, where the sequence of s's is aligned
% with the appropriate k's, according to which segment contains the
% x coordinate.  Note that
% s := (x - r_i)/(r_{i+1} - r_i).
%
% Version:        Changes:
% --------        -------------
% 05.12.2016      Original code.
% 04.07.2017      Adapted for complex step.
% 08.08.2017      Added dS/dx as an output.
%
% Version:        Verification:
% --------        -------------
% 05.12.2016      Memo AN 16.12.70
% 04.07.2017      
% 08.08.2017      
%
% Inputs:
% -------
% rhat            : A list of x coordinates of the control points.
% x               : A list of x coordinates of the points.
%
% Outputs:
% --------
% S               : A matrix relating z(x) = S(x)*k.
% dSdx            : dz/dx = dS/dx*k.
%

Nr = size(rhat,1);
Nx = size(x,1);

% I can get just what I need from the interp1 function.  Let the X
% values be the rhat's, and the Y's be the indices.  The outputs
% are then floating point numbers, the integer part of which 
% indicates the index of the interval (position in the S matrix),
% and the decimal part is precisely the value of s that I'm looking
% for.
y = interp1(real(rhat),[1:Nr].',real(x),'extrap');

inter = floor(y);
inter(inter<1) = 1;            % Prevent going outside the
inter(inter>Nr-1) = Nr - 1;    % valid reference interval.
ind = 4*(inter-1);
den = rhat(inter+1) - rhat(inter);
ss = (x - rhat(inter))./den;

ii = [[1:Nx] [1:Nx] [1:Nx] [1:Nx]].';
jj = [ind.'+1 ind.'+2 ind.'+3 ind.'+4].';

qq = [ones(1,Nx) ss.' (ss.').^2 (ss.').^3];

S = sparse(ii,jj,qq,Nx,4*(Nr-1));

qq = [zeros(1,Nx)./(den.') ...
       ones(1,Nx)./(den.') ...
         (2*ss.')./(den.') ...
    (3*(ss.').^2)./(den.')];

dSdx = sparse(ii,jj,qq,Nx,4*(Nr-1));