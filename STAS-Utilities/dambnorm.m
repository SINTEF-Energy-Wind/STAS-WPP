function d = dambnorm (a,b)
%           d   a - b
% Compute   -- -------
%           da |a - b|

N = size(a,1);
d = zeros(N,N);

amb = a - b;
namb = sqrt((amb.')*amb);
namb2 = namb^2;

for jj = 1:N

   d(:,jj) = -amb*amb(jj);
   d(jj,jj) = d(jj,jj) + namb^2;
   d(:,jj) = d(:,jj)/(namb^3);

end
