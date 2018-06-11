function bsh = bladeModeShape (s,ret,shape)
% 
% Extract the blade elastic mode shapes from the mode shape matrix.
%
% Version:        Changes:
% --------        -------------
% 16.04.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 16.04.2018      Checked output from a sample calculation.
%

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

Nret = size(ret,1);
ind = [1:Nret].';

imod = 6 + s.foundation.Nmod + s.tower.Nmod ...
     + s.nacelle.Nmod + s.driveshaft.Nmod;
Nmod = s.blade(1).Nmod;
DOF1 = ind(ret==idofs(6)+7);
DOF2 = ind(ret==idofs(7));

% Extract the flapwise component at the element location.
sh1 = [sparse(1,Nmod);shape([DOF1+2:6:DOF2-9],imod+[1:Nmod])];
sh2 = shape([DOF1+2:6:DOF2-3],imod+[1:Nmod]);

bsh = 0.5*(sh1 + sh2);

