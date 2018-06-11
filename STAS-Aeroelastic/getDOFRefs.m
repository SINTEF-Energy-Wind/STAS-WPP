function [idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s)
%
% Version:        Changes:
% --------        -------------
% 26.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 26.10.2017      Frequently used in the code.  No errors observed.
%

idofs = zeros(8,1);
idofm = zeros(8,1);
inods = zeros(8,1);
inodm = zeros(8,1);

% Global to foundation.
idofm(1) = 0;
inodm(1) = 0;
idofs(1) = 0;
inods(1) = 1;

% Foundation to tower.
idofm(2) = s.foundation.Ndof - 6;
inodm(2) = s.foundation.Nnod;
idofs(2) = idofm(2) + 6;
inods(2) = inodm(2) + 1;

% Tower to nacelle.
idofm(3) = idofm(2) + s.tower.Ndof;
inodm(3) = inodm(2) + s.tower.Nnod;
idofs(3) = idofm(3) + 6;
inods(3) = inodm(3) + 1;

% Nacelle to driveshaft, rear bearing.
idofm(4) = idofm(3) + 6*(s.Nev + s.Ner + 1);
inodm(4) = inodm(3) + s.Nev + s.Ner + 1;
idofs(4) = idofs(3) + s.nacelle.Ndof;
inods(4) = inods(3) + s.nacelle.Nnod;

% Nacelle to driveshaft, front bearing.
idofm(5) = idofm(3) + s.nacelle.Ndof;
inodm(5) = inodm(3) + s.nacelle.Nnod;
idofs(5) = idofs(4) + 6*s.Ned;
inods(5) = inods(4) + s.Ned;

% Driveshaft to blade 1.
idofm(6) = idofs(4) + s.driveshaft.Ndof - 18;
inodm(6) = inods(4) + s.driveshaft.Nnod - 3;
idofs(6) = idofs(4) + s.driveshaft.Ndof;
inods(6) = inods(4) + s.driveshaft.Nnod;

% Driveshaft to blade 2.
idofm(7) = idofm(6) + 6;
inodm(7) = inodm(6) + 1;
idofs(7) = idofs(6) + s.blade(1).Ndof;
inods(7) = inods(6) + s.blade(1).Nnod;

% Driveshaft to blade 3.
idofm(8) = idofm(7) + 6;
inodm(8) = inodm(7) + 1;
idofs(8) = idofs(7) + s.blade(2).Ndof;
inods(8) = inods(7) + s.blade(2).Nnod;

Ndof = idofs(8) + s.blade(3).Ndof;