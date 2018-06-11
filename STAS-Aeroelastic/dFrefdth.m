function dF = dFrefdth (TTh,TBB0,dTBB0,d2TBB0,TB0g)
%
% Version:        Changes:
% --------        -------------
% 26.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 26.10.2017      Verified by complex step.
%

dF = zeros(3,3*3*3);

%d2TBB0 = d2Tdth2 (TTh,TBB0,dTBB0);

for idof = 1:3

   ic9 = 9*(idof-1);
   ic3 = 3*(idof-1);

   for jj = 1:3

      jc3 = 3*(jj-1);

      dF(:,ic9+jc3+[1:3]) = TB0g*d2TBB0(:,ic9+jc3+[1:3])*(TB0g*TBB0).' ...
                          + TB0g*dTBB0(:,jc3+[1:3])*(TB0g*dTBB0(:,ic3+[1:3])).';
      
   end

end