function [b1,b2,b3] = MBCindices_Nmd (Nmd,imdofs)

Nmb = imdofs(7) - imdofs(6);

ind = [1:Nmd].';
for ib = 1:3
   ii = zeros(Nmd,1);
   ii(imdofs(5+ib)+[1:Nmb]) = 1;
   ii = logical(ii);
   if (ib == 1)
      b1 = ind(ii);
   elseif (ib == 2)
      b2 = ind(ii);
   elseif (ib == 3)
      b3 = ind(ii);
   end
end

