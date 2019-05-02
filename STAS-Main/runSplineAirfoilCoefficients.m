% Generate splined airfoil coefficients.

nm = 'Tjaereborg'

clear;

% First column is aoa (deg), then trios of Cl,Cd,Cm for each airfoil type.
load(['ClCdCm_' nm '.txt']);
eval(['ClCdCm = ClCdCm_' nm]);

aoas = ClCdCm(:,1)*pi/180;

[A,B] = splineCoefficients (aoas);
k = (A\B)*ClCdCm(:,2:size(ClCdCm,2));
kd = zeros(size(k));
Nk = size(k,1);
for jj = 1:Nk/4
   irow = 4*(jj-1);
   kd(irow+[1:3],:) = k(irow+[2:4],:).*[1;2;3];
end
kfoils = [k kd];

% Find the zero-lift angles-of-attack, needed for dynamic stall.
% These values may need manual modification in the event that
% the airfoils do not behave as normal airfoils.  This can be
% the case for very thick inboard profiles.
Nfoils = (size(ClCdCm,2) - 1)/3;
aoazs = zeros(Nfoils,1);
for ifoil = 1:Nfoils

   wt        = zeros(Nfoils,1);
   wt(ifoil) = 1;

   [aoazs(ifoil),fval,exFlg,outp] = ...
      fzero (@(aoa) getLift (aoas,kfoils,wt,aoa),[-pi/4;pi/4]);

end

% Valid for Gnu Octave only.
save('-binary',['aoas_' nm '.bin'],'aoas');
save('-binary',['kfoils_' nm '.bin'],'kfoils');
save('-binary',['aoazs_' nm '.bin'],'aoazs');

%weights = [0.5 0.5 0 0 0 0].';
%aoa = [-180:1:180]'*(pi/180);
%C = airfoilCoefficientsSpline (aoas,kfoils,weights,aoa);
%[aoa*180/pi C]
