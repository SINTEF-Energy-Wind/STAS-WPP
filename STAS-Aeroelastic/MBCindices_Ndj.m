function [b1,b2,b3] = MBCindices_Ndj (Ndj,idofs)

Ndb = idofs(7) - idofs(6);

ind = [1:Ndj].';
for ib = 1:3
   ii = zeros(Ndj,1);
   ii(idofs(5+ib)+[1:Ndb]) = 1;
   ii = logical(ii);
   if (ib == 1)
      b1 = ind(ii);
   elseif (ib == 2)
      b2 = ind(ii);
   elseif (ib == 3)
      b3 = ind(ii);
   end
end



 