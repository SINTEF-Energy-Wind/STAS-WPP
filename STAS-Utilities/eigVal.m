function [slap,shp,ifrq] = eigVal (A)
%
% Compute and print eigenvalues and eigenvectors of a matrix.
%

[Diag,Perm,Abal] = balance(A,'noperm');
Neig = size(Abal,1);
DD = sparse([1:Neig],[1:Neig],Diag,Neig,Neig);
iDD = inv(DD);
[shpb,lam] = eig(Abal);
shp = DD*shpb*iDD;
f1 = diag(lam);
[si,ifrq] = sort(imag(f1));
slap = f1(ifrq);



zeta = zeros(Neig,1);
for j = 1:Neig
   if (abs(imag(slap(j))) > 1e-6*abs(real(slap(j))))
      R = real(slap(j))/imag(slap(j));
   else
      R = real(slap(j)); % 0;
   end
   zeta(j) = sqrt((R^2)/(1 + R^2));
   if (real(slap(j)) > 0)
      zeta(j) = -zeta(j);
   end

   % Scale the mode shapes to unit norm length.  (Note that the
   % shapes are not ordered the same as the sorted slap... but
   % it doesn't matter, we're looping through all the shapes.)
   shp(:,j) = shp(:,j)/norm(shp(:,j),2);

end

[[1:Neig]' zeta real(slap)./(2*pi) imag(slap)./(2*pi)]


