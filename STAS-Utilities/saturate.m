function [y,Du,D2u] = saturate (u,params)
%
% A unit saturation function, with splined transitions.  The saturation region
% is given a slight slope for reasons of numerical stability when integrated
% with the larger code.
%
% Version:        Changes:
% --------        -------------
% 11.01.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 11.01.2019      Inspected outputs.
%
% Inputs:
% -------
% u               : x-axis values.
% params          : xa,xb, the u intersections at the beginning and end of
%                   the transition curve.
%
% Outputs:
% --------
% y, Du = dy/du, D2u = d^2y/du^2.

Nu = size(u,1);
y = zeros(Nu,1);
Dvec = zeros(Nu,1);
D2vec = zeros(Nu,1);

xa = params(1);
xb = params(2);

xa2 = xa^2;
xa3 = xa2*xa;
xb2 = xb^2;
xb3 = xb2*xb;

smin = 0.0;  % 0.05;

ind = (real(u) <= -xb);
y(ind) = -1 + smin*(u(ind) + xb);
Dvec(ind) = smin;

ind = (real(u) >= xb);
y(ind) = 1 + smin*(u(ind) - xb);
Dvec(ind) = smin;

ind = (real(u) >= -xa) & (real(u) <= xa);
y(ind) = u(ind);
Dvec(ind) = 1;

ind = (real(u) > xa) & (real(u) < xb);
M = [1   xa   xa2   xa3; ...
     0   1   2*xa 3*xa2; ...
     1   xb   xb2   xb3; ...
     0   1   2*xb 3*xb2];
b = [xa;1;1;smin];             % Give it a slight slope for numerical stability.
c = M\b;
y(ind) = c(1) + u(ind).*(c(2) + u(ind).*(c(3) + u(ind)*c(4)));
Dvec(ind) = c(2) + u(ind).*(2*c(3) + u(ind)*3*c(4));
D2vec(ind) = 2*c(3) + u(ind)*6*c(4);

ind = (real(u) < -xa) & (real(u) > -xb);
y(ind) = -c(1) + u(ind).*(c(2) + u(ind).*(-c(3) + u(ind)*c(4)));
Dvec(ind) = c(2) + u(ind).*(-2*c(3) + u(ind)*3*c(4));
D2vec(ind) = -2*c(3) + u(ind)*6*c(4);

ii = [1:Nu];
jj = [1:Nu];
Du = sparse(ii,jj,Dvec,Nu,Nu);
D2u = sparse(ii,jj,D2vec,Nu,Nu);



