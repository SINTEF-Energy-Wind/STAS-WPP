function d2 = d2axbdab (a,b)
%             d^2    a x b
% Compute   ------- -------
%           dai dbj |a x b|
%
% 08.11.2017: verified by finite difference.

N = 3;
d2 = zeros(N,N*N);

axb = cross(a,b);
axbTaxb = axb.'*axb;

for ii = 1:N

   icN2 = N*(ii-1);

   ei = zeros(N,1);
   ei(ii) = 1;
   eixb = cross(ei,b);

   for jj = 1:N

      ej = zeros(N,1);
      ej(jj) = 1;
      eixej = cross(ei,ej);
      axej = cross(a,ej);

      d2(:,icN2+jj) = -3*(eixb*axbTaxb - axb*((eixb.')*axb))     ...
                    *    ((axej.')*axb)/(axbTaxb^(5/2))          ...
                    + (eixej*axbTaxb + 2*eixb*((axej.')*axb)     ...
                    -  axej*((eixb.')*axb) - axb*((eixej.')*axb) ...
                    -  axb*((eixb.')*axej))/(axbTaxb^(3/2));

   end

end

