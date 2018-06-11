function d2 = d2ambnorm (a,b)
%             d^2    a - b
% Compute   ------- -------
%           dai daj |a - b|
%
% 08.11.2017: verified by complex step.

N = size(a,1);
d2 = zeros(N,N*N);

amb = a - b;
ambTamb = (amb.')*amb;

for ii = 1:N

   icN2 = N*(ii-1);

   ei = zeros(N,1);
   ei(ii) = 1;

   for jj = 1:N

      ej = zeros(N,1);
      ej(jj) = 1;

      d2(:,icN2+jj) = -3*(ei*ambTamb - amb*amb(ii))*amb(jj)/(ambTamb^(5/2)) ...
                    + (2*ei*amb(jj) - ej*amb(ii) - amb*((ei.')*ej))/(ambTamb^(3/2));

   end

end


