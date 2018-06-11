function [TmT,dTmT] = dTmeT (mes,TsB,dTsB)
%
% A reduced version of d2TmeT for when the second derivative is
% not needed, in particular building the nonlinear equations.
%
% Version:        Changes:
% --------        -------------
% 23.11.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 23.11.2017      Same result as d2TmeT.
%

TsBe = [TsB zeros(3,9);            ...
        zeros(3,3) TsB zeros(3,6); ...
        zeros(3,6) TsB zeros(3,3); ...
        zeros(3,9) TsB];
TBse = TsBe.';

TmT = TsBe*mes*TBse;
dTmT = zeros(12,12*18);
for ii = 1:12

   ic12 = 12*(ii+6-1);  % First six groups are zero.
   ic3 = 3*(ii-1);

   dTi = [dTsB(:,ic3+[1:3]) zeros(3,9);            ...
          zeros(3,3) dTsB(:,ic3+[1:3]) zeros(3,6); ...
          zeros(3,6) dTsB(:,ic3+[1:3]) zeros(3,3); ...
          zeros(3,9) dTsB(:,ic3+[1:3])];
   mat = dTi*mes*TBse;
   dTmT(:,ic12+[1:12]) = mat + mat.';

end
dTmT = sparse(dTmT);

