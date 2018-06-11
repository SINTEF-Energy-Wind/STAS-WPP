function d2 = d2axbnorm (a,b)
%             d^2    a x b
% Compute   ------- -------
%           dai daj |a x b|
%
% 08.11.2017: verified by complex step.

N = 3;
d = zeros(N,N*N);

axb = cross(a,b);
axbTaxb = (axb.')*axb;

for ii = 1:N

   icN2 = N*(ii-1);

   ei = zeros(N,1);
   ei(ii) = 1;
   eixb = cross(ei,b);

   for jj = 1:N

      ej = zeros(N,1);
      ej(jj) = 1;
      ejxb = cross(ej,b);

      d2(:,icN2+jj) = -3*(eixb*axbTaxb - axb*((eixb.')*axb))*((ejxb.')*axb) ...
                    /    (axbTaxb^(5/2))                                    ...
                    + (2*eixb*((ejxb.')*axb) - ejxb*((eixb.')*axb)          ...
                    -  axb*((eixb.')*ejxb))                                 ...
                    / (axbTaxb^(3/2));

   end

end
