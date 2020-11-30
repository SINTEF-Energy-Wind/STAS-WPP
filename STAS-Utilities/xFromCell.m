function x = xFromCell (ic,xb)
%
% Return the cell states for a given cell ID.
%
% Version:        Changes:
% --------        -------------
% 17.10.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 17.10.2018      
%
% Inputs:
% -------
% ic              : nearest cell id.
% xb              : N-by-3 array containing xmin, dx, xmax for each x;
%                   that is, the possible values.
%
% Outputs:
% --------
% x               : the vector of N states.

icr = ic;
Ndim = size(xb,1);
Nvec = size(ic,1);
Nx = round((xb(:,3) - xb(:,1))./xb(:,2)) + 1;
ix = zeros(Ndim,Nvec);
for idim = 1:Ndim

   if (idim == Ndim)
      ix(idim,:) = icr;
   else
      pr = prod(Nx(idim+1:Ndim));
      n = floor(icr/pr);
      rem = mod(icr,pr);
      n(rem==0) = n(rem==0) - 1;
      ix(idim,:) = n + 1;
%[idim pr icr(1) floor(icr(1)/pr) n ix(idim,1)]
      icr = icr - pr*n;
   end

end
%ix
x = xb(:,1) + (ix - 1).*xb(:,2);
