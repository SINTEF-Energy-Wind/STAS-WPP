function [d,ind] = cluster (x,nn)
%
% Find the nearest neighbors.
%
% Version:        Changes:
% --------        -------------
% 30.06.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 30.06.2020      Tested on simple cases.
%
% Inputs:
% -------
% x               : Positions, Ndim-by-Np.
% nn              : Number of neighbors to return.
%
% Outputs:
% --------
% d               : Matrix of distances, nn-by-Np.

% This is an inefficient algorithm when it comes to storage and
% operations, but it is efficient to program, and it is vectorized
% in a way that is intended for small dimensions but potentially
% large numbers of points.

Ndim = size(x,1);
Np = size(x,2);

dx = zeros(Np,Np);
for idim = 1:Ndim
   dx = dx + (x(idim,:).' - x(idim,:)).^2;
end
dx = sqrt(dx);
[s,k] = sort(dx);
d = s(1:nn,:);
ind = k(1:nn,:);
