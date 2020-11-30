function [ic,wc] = nearestCells (x,xb)
%
% Return the nearest cell indices, for the cells surrounding each
% entry in a vector of states.  Also the weights for a linear
% interpolation.
%
% Version:        Changes:
% --------        -------------
% 21.05.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 21.05.2019      
%
% Inputs:
% -------
% x               : Ndim-by-Nvecs array containing the vectors.
% xb              : Ndim-by-3 array containing xmin, dx, xmax for each dimension;
%                   that is, the possible values.
%
% Outputs:
% --------
% ic              : nearest cell ids, 2^Ndim-by-Nvecs.
% wc              : weights associated with each of the ic cells.

Ndim = size(x,1);
Nvecs = size(x,2);
Nx = round((xb(:,3) - xb(:,1))./xb(:,2)) + 1;
vx  = (x - xb(:,1))./xb(:,2);
ixx = min(max(vx + 1,1),Nx);
ixL = min(max(floor(vx) + 1,1),Nx-1);
wL  = 1 - (ixx - ixL);

ic = zeros(2^Ndim,Nvecs);
wc = zeros(2^Ndim,Nvecs);

for jj = 1:2^Ndim

   b = dec2bin(jj-1,Ndim);
   ix = ixL;

   ind = logical(bin2dec(b(:)));
   ix(ind,:) = ix(ind,:) + 1;  % The values that are 1 in binary, make "high".

   wx = wL;
   wx(ind,:) = 1 - wL(ind,:);
   wc(jj,:) = prod(wx,1);

   for idim = 1:Ndim

      if (idim == Ndim)
         ic(jj,:) = ic(jj,:) + ix(idim,:);
      else
         ic(jj,:) = ic(jj,:) + prod(Nx(idim+1:Ndim))*(ix(idim,:)-1);
      end      

   end
end

ic = min(max(ic,1),prod(Nx));
