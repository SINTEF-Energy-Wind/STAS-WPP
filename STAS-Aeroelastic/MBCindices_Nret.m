function [b1,b2,b3] = MBCindices_Nret (Ndj,idofs,ret,slv)

Nret = size(ret,1);
Ndb = idofs(7) - idofs(6);

indr = [1:Nret].';
for ib = 1:3
   ii = zeros(Ndj,1);
   ii(idofs(5+ib)+[1:Ndb]) = 1;
   [ip,rr,cr] = partitionMatrix(ii,slv,[]);
   ir = ip(1:Nret);
   ir = logical(ir);
   if (ib == 1)
      b1 = indr(ir);
   elseif (ib == 2)
      b2 = indr(ir);
   elseif (ib == 3)
      b3 = indr(ir);
   end
end

