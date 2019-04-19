function [TmT,dTmT,d2TmT] = d2TmeT (mes,TsB,dTsB,d2TsB)
%
% Version:        Changes:
% --------        -------------
% 13.11.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 13.11.2017      Verified using complex step derivatives.
%

TsBe = [TsB zeros(3,9);            ...
        zeros(3,3) TsB zeros(3,6); ...
        zeros(3,6) TsB zeros(3,3); ...
        zeros(3,9) TsB];
TBse = TsBe.';

TmT = TsBe*mes*TBse;
dTmT = zeros(12,12*18);
d2TmT = zeros(12,12*18*18);
mT = mes*TBse;
for ii = 1:12

   ic216 = 216*(ii+6-1); % First six groups are zero.
   ic36 = 36*(ii-1);
   ic12 = 12*(ii+6-1);
   ic3 = 3*(ii-1);

   dTi = [dTsB(:,ic3+[1:3]) zeros(3,9);            ...
          zeros(3,3) dTsB(:,ic3+[1:3]) zeros(3,6); ...
          zeros(3,6) dTsB(:,ic3+[1:3]) zeros(3,3); ...
          zeros(3,9) dTsB(:,ic3+[1:3])];
   dTm = dTi*mes;
   mat = dTm*TBse;
   dTmT(:,ic12+[1:12]) = mat + mat.';

   for jj = 1:12

      jc12 = 12*(jj+6-1);      
      jc3 = 3*(jj-1);

      dTj = [dTsB(:,jc3+[1:3]) zeros(3,9);            ...
             zeros(3,3) dTsB(:,jc3+[1:3]) zeros(3,6); ...
             zeros(3,6) dTsB(:,jc3+[1:3]) zeros(3,3); ...
             zeros(3,9) dTsB(:,jc3+[1:3])];
      d2T = [d2TsB(:,ic36+jc3+[1:3]) zeros(3,9);            ...
             zeros(3,3) d2TsB(:,ic36+jc3+[1:3]) zeros(3,6); ...
             zeros(3,6) d2TsB(:,ic36+jc3+[1:3]) zeros(3,3); ...
             zeros(3,9) d2TsB(:,ic36+jc3+[1:3])];
      mat = d2T*mT + dTm*(dTj.');
      d2TmT(:,ic216+jc12+[1:12]) = mat + mat.';

   end

end
dTmT = sparse(dTmT);
d2TmT = sparse(d2TmT);
