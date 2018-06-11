function d2F = d2Fdth2 (TTh,th,TB0g,Tn0B)
%
% Version:        Changes:
% --------        -------------
% 06.11.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 06.11.2017      Verified by finite difference and complex step.
%

d2F = zeros(3,3*3*6*6);  % Initial derivative wrt th, higher derivs wrt Phi,th.

[Tnn0,dTnn0] = dTdth (th);
d2Tnn0 = d2Tdth2 (th,Tnn0,dTnn0);
d3Tnn0 = d3Tdth3 (th,Tnn0,dTnn0,d2Tnn0);
[TBB0,dTBB0] = dTdth (TTh);
d2TBB0 = d2Tdth2 (TTh,TBB0,dTBB0);

% d2Fi/(dPhij dPhik)
for jj = 1:3

   jc54 = 54*(jj-1);
   jc9  =  9*(jj-1);
   jc3  =  3*(jj-1);

   for kk = 1:3

      kc9 = 9*(kk-1);
      kc3 = 3*(kk-1);

      for ii = 1:3

         ic3 = 3*(ii-1);

         term1 = d2TBB0(:,jc9+kc3+[1:3])*Tn0B*dTnn0(:,ic3+[1:3]) ...
               * ((TBB0*Tn0B*Tnn0).');
         term2 = dTBB0(:,jc3+[1:3])*Tn0B*dTnn0(:,ic3+[1:3]) ...
               * ((dTBB0(:,kc3+[1:3])*Tn0B*Tnn0).');
         term3 = dTBB0(:,kc3+[1:3])*Tn0B*dTnn0(:,ic3+[1:3]) ...
               * ((dTBB0(:,jc3+[1:3])*Tn0B*Tnn0).');
         term4 = TBB0*Tn0B*dTnn0(:,ic3+[1:3]) ...
               * ((d2TBB0(:,jc9+kc3+[1:3])*Tn0B*Tnn0).');

         d2F(:,jc54+kc9+ic3+[1:3]) = ...
               TB0g*(term1 + term2 + term3 + term4)*(TB0g.');

      end

   end

end

% d2Fi/(dPhij dthk)
for jj = 1:3

   jc54 = 54*(jj-1);
   jc9  =  9*(jj-1);
   jc3  =  3*(jj-1);

   for kk = 1:3

      kc9a = 9*(kk+3-1);
      kc9 = 9*(kk-1);
      kc3 = 3*(kk-1);

      for ii = 1:3

         ic3 = 3*(ii-1);

         term1 = dTBB0(:,jc3+[1:3])*Tn0B*d2Tnn0(:,kc9+ic3+[1:3]) ...
               * ((TBB0*Tn0B*Tnn0).');
         term2 = dTBB0(:,jc3+[1:3])*Tn0B*dTnn0(:,ic3+[1:3]) ...
               * ((TBB0*Tn0B*dTnn0(:,kc3+[1:3])).');
         term3 = TBB0*Tn0B*d2Tnn0(:,kc9+ic3+[1:3]) ...
               * ((dTBB0(:,jc3+[1:3])*Tn0B*Tnn0).');
         term4 = TBB0*Tn0B*dTnn0(:,ic3+[1:3]) ...
               * ((dTBB0(:,jc3+[1:3])*Tn0B*dTnn0(:,kc3+[1:3])).');

         d2F(:,jc54+kc9a+ic3+[1:3]) = ...
               TB0g*(term1 + term2 + term3 + term4)*(TB0g.');

      end

   end

end

% d2Fi/(dthj dPhik)
for jj = 1:3

   jc54a = 54*(jj+3-1);
   jc9a  =  9*(jj+3-1);

   for kk = 1:3

      kc9  =  9*(kk-1);
      kc54 = 54*(kk-1);

      for ii = 1:3

         ic3 = 3*(ii-1);

         d2F(:,jc54a+kc9+ic3+[1:3]) = d2F(:,kc54+jc9a+ic3+[1:3]);

      end

   end

end


% d2Fi/(dthj dthk)
for jj = 1:3

   jc54a = 54*(jj+3-1);
   jc27  = 27*(jj-1);
   jc9   =  9*(jj-1);
   jc3   =  3*(jj-1);

   for kk = 1:3

      kc9a = 9*(kk+3-1);
      kc9 = 9*(kk-1);
      kc3 = 3*(kk-1);

      for ii = 1:3

         ic3 = 3*(ii-1);

         term1 = TBB0*Tn0B*d3Tnn0(:,jc27+kc9+ic3+[1:3]) ...
               * ((TBB0*Tn0B*Tnn0).');
         term2 = TBB0*Tn0B*d2Tnn0(:,jc9+ic3+[1:3]) ...
               * ((TBB0*Tn0B*dTnn0(:,kc3+[1:3])).');
         term3 = TBB0*Tn0B*d2Tnn0(:,kc9+ic3+[1:3]) ...
               * ((TBB0*Tn0B*dTnn0(:,jc3+[1:3])).');
         term4 = TBB0*Tn0B*dTnn0(:,ic3+[1:3]) ...
               * ((TBB0*Tn0B*d2Tnn0(:,jc9+kc3+[1:3])).');

         d2F(:,jc54a+kc9a+ic3+[1:3]) = ...
               TB0g*(term1 + term2 + term3 + term4)*(TB0g.');

      end

   end

end
