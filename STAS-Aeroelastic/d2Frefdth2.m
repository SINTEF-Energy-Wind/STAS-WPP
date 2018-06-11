function d2F = d2Frefdth2 (TTh,TBB0,dTBB0,d2TBB0,d3TBB0,TB0g)
%
% Version:        Changes:
% --------        -------------
% 06.11.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 06.11.2017      Verified by finite difference using dFrefdth.m.
%

d2F = zeros(3,3*3*3*3);

%d3TBB0 = d3Tdth3 (TTh,TBB0,dTBB0,d2TBB0);

for jj = 1:3

   jc27 = 27*(jj-1);
   jc9  =  9*(jj-1);
   jc3  =  3*(jj-1);

   for kk = 1:3

      kc9 = 9*(kk-1);
      kc3 = 3*(kk-1);

      for ii = 1:3

         ic3 = 3*(ii-1);

         d2F(:,jc27+kc9+ic3+[1:3]) =                                             ...
                  TB0g*(d3TBB0(:,jc27+kc9+ic3+[1:3])*( TBB0.')                   ...
                +       d2TBB0(:,jc9+ic3+[1:3])     *( dTBB0(:,kc3+[1:3]).')     ...
                +       d2TBB0(:,kc9+ic3+[1:3])     *( dTBB0(:,jc3+[1:3]).')     ...
                +        dTBB0(:,ic3+[1:3])         *(d2TBB0(:,jc9+kc3+[1:3]).') ...
                       )*(TB0g.');

      end
      
   end

end