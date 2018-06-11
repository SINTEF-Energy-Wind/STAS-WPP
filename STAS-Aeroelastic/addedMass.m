function Madd = addedMass (Nnf,Nwater,wel,Lel,rhow,CmA,Qh)
%
% Hydrodynamic added mass.  The force on the circular foundation is
% F = - rhow (Cm A) L d^2x/dt^2 = - M' d^2x/dt^2 = -M' Qh' d^2q/dt^2.
% Applying this force at the nodes requires the operation G = Qh F,
% so the result, after moving to the LHS of the equations, is
% Madd = Qh M' Qh'.
%
% Version:        Changes:
% --------        -------------
% 27.11.2017      Code adapted from earlier addedMass.m.  Fixed the
%                 omission of Qhat and changed inputs.
%
% Version:        Verification:
% --------        -------------
% 27.11.2017      
%

wam = wel(1:Nwater);  % Only full-time submerged elements get added mass.
dvec = zeros(6*Nnf,1);
idof = [6*(wam-1)+1;6*(wam-1)+2];
iel = [wam;wam];
C = [CmA(1:Nwater);CmA(1:Nwater)];
dvec(idof) = 0.5*rhow*C.*Lel(iel);  % Half-length.
idof = [6*wam+1;6*wam+2];
dvec(idof) = dvec(idof) + 0.5*rhow*C.*Lel(iel);
idof = [1:6*Nnf].';
Mp = sparse(idof,idof,dvec,6*Nnf,6*Nnf);
Madd = Qh*Mp*(Qh.');


