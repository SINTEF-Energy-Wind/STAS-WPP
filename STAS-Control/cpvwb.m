function [cp,ct,dcp,dct] = cpvwb (cpct,R,Vinf,W,bet)
%
% Version:        Changes:
% --------        -------------
% 08.02.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 08.02.2019      Checked sample cases.  Derivatives verified by complex step.
%

del = eps^0.5;
lam = R*W/maxc(Vinf,del);

% Find the indices through 2D interpolation.
Nl = cpct(1,1);
Nb = cpct(1,2);
lams = cpct(1+[1:Nb:Nb*(Nl-1)+1],1);
bets = cpct(2:Nb+1,2);

ilam = floor (interp1 (lams,[1:Nl],real(lam),'extrap'));
ilam = max(min(ilam,Nl-1),1);

ibet = floor (interp1 (bets,[1:Nb],real(bet),'extrap'));
ibet = max(min(ibet,Nb-1),1);

ind = 1 + [(ilam-1)*Nb+ibet; ilam*Nb+ibet; ilam*Nb+ibet+1; (ilam-1)*Nb+ibet+1];

xv  = [cpct(ind,1) cpct(ind,2)];

y   = cpct(ind,3);
y1  = cpct(ind,4);
y2  = cpct(ind,5);
y12 = cpct(ind,6);
[cp,dcp1] = bicub ([lam bet],y,y1,y2,y12,xv);

y   = cpct(ind,7);
y1  = cpct(ind,8);
y2  = cpct(ind,9);
y12 = cpct(ind,10);
[ct,dct1] = bicub ([lam bet],y,y1,y2,y12,xv);

dcp = [dcp1(1)*(-lam/maxc(Vinf,del)); dcp1(1)*(R/maxc(Vinf,del)); dcp1(2)];
dct = [dct1(1)*(-lam/maxc(Vinf,del)); dct1(1)*(R/maxc(Vinf,del)); dct1(2)];
