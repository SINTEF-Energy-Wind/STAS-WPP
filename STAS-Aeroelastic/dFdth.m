function dF = dFdth (Tnn0,dTnn0,d2Tnn0,TBB0,dTBB0,TB0g,Tn0B)
%
% Version:        Changes:
% --------        -------------
% 26.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 26.10.2017      Verified by complex step and finite difference.
%

dF = zeros(3,3*3*6);

Tn0g = TB0g*TBB0*Tn0B;

for idof = 1:3   % First the three body rotations.

   ic9 = 9*(idof-1);
   ic3 = 3*(idof-1);

   TdT = TB0g*dTBB0(:,ic3+[1:3])*Tn0B;

   for jj = 1:3
      jc3 = 3*(jj-1);
      dF(:,ic9+jc3+[1:3]) = TdT*dTnn0(:,jc3+[1:3])*(Tn0g*Tnn0).' ...
                          + Tn0g*dTnn0(:,jc3+[1:3])*(TdT*Tnn0).';
   end

end

for idof = 1:3   % Then the nodal rotations.

   jc9 = 9*(idof+3-1);
   ic9 = 9*(idof-1);
   ic3 = 3*(idof-1);

   for jj = 1:3
      jc3 = 3*(jj-1);
      dF(:,jc9+jc3+[1:3]) = Tn0g*d2Tnn0(:,ic9+jc3+[1:3])*(Tn0g*Tnn0).' ...
                          + Tn0g*dTnn0(:,jc3+[1:3])*(Tn0g*dTnn0(:,ic3+[1:3])).';
   end

end