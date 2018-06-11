function [shape,freq] = modes (Nmod,M,K)
%
% Version:        Changes:
% --------        -------------
% 04.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 04.12.2017      Used in bodyModes, same verification status.
%                 Experience indicates that this code produces the 
%                 correct results.
%

shape = zeros(size(M,1),Nmod);
freq  = zeros(Nmod,1);

[shp,lam] = eig(K,M);
f1 = diag(lam);
[w2,ifreq] = sort(f1);
shp = shp(:,ifreq(1:Nmod));
w2 = w2(1:Nmod);

%[shp,lam] = eigs(K,M,Nmod,'sm');
%w2 = diag(lam);

for ieig = 1:Nmod
   % "max" here is OK for complex step, as it is used only to get
   % the index imax.
   %
   % Arrange things such that the maximum magnitude entry in the
   % normalized shape is positive.
   [scale,imax] = max(real(absc(shp(:,ieig))));
   smax = sign(real(shp(imax,ieig)));
%   shape(:,ieig) = shp(:,ieig)/shp(imax,ieig);
   shape(:,ieig) = smax*shp(:,ieig)/sqrt((shp(:,ieig).')*shp(:,ieig));
   freq(ieig) = sqrt(w2(ieig))/(2*pi);

%f = sqrt(w2(ieig))/(2*pi)
%printVec(shape(:,ieig))

end

