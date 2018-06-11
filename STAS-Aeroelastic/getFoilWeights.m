function fwt = getFoilWeights (tcfoils,tca)
% Compute foil weights from t/c.  This works with complex step.
%
% Version:        Changes:
% --------        -------------
% 01.11.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 01.11.2017      Checked results for various analyses including IEA
%                 aero optimization studies.
%

Nfoils = size(tcfoils,1);
Neb = size(tca,1);

y = interp1(real(tcfoils),[1:Nfoils],real(tca),'extrap');
inter = floor(y);
inter(inter<1) = 1;
inter(inter>Nfoils-1) = Nfoils-1;

den = tcfoils(inter+1) - tcfoils(inter);
wt = (tca - tcfoils(inter))./den;

ii = [inter;inter+1];
jj = [[1:Neb] [1:Neb]].';
ss = [(1-wt);wt];
fwt = sparse(ii,jj,ss,Nfoils,Neb);
