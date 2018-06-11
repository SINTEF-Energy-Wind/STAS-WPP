function d = daxbnorm (a,b)
%           d   a x b
% Compute   -- -------
%           da |a x b|
%
% 08.11.2017: verified by finite difference.

N = 3;
d = zeros(N,N);

axb = cross(a,b);
naxb = sqrt((axb.')*axb);
naxb3 = naxb^3;

for jj = 1:N

   ee = zeros(N,1);
   ee(jj) = 1;
   exb = cross(ee,b);
   d(:,jj) = exb/naxb - axb*((exb.')*axb)/naxb3;

end
