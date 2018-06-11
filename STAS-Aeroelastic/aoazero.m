function aoaz = aoazero (aoas,kfoils,foilwt)
%
% Find the zero-lift angle-of-attack for an airfoil.
%
% Version:        Changes:
% --------        -------------
% 14.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 14.12.2017      Checked that this works for the DTU 10 MW airfoils.
%                 Does not work for the inboard "airfoils" that have
%                 a highly unusual lift-aoa profile.
%

Nel = size(foilwt,2);

aoaz = zeros(Nel,1);
for iel = 1:Nel

   [aoaz(iel),fval,stat,outp] =                                     ...
               fzero (@(aoa) getCL (aoas,kfoils,foilwt(:,iel),aoa), ...
                      [-pi/4;pi/4]);

end