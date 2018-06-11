function Psi = aeroPsi (a,rp,bsh)
%
% Compute the aerodynamic mode shapes or reduced basis.
%
% Version:        Changes:
% --------        -------------
% 25.04.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 25.04.2018      Simple function, checked output for sample cases.
%
% Inputs:
% -------
% a               : Aero data structure.
% bsh             : For a.icp < 0, blade mode shapes.
% rp              : For a.icp > 0, projected radii of blade elements.
%
% Outputs:
% --------
% Psi             : Nxa-by-Nra matrix.

Nae = a.Neb*a.Nb;
Nxa = 7*Nae;
Neb = a.Neb;

if (a.icp(1) == 0)

   % No aero mode reduction.
   Psi = speye(Nxa);

else

   Ncp = size(a.icp,1);

   i7a = [1:7:7*Neb-6].';
   i7b = [2:7:7*Neb-5].';
   i7c = [3:7:7*Neb-4].';
   i7d = [4:7:7*Neb-3].';
   i7e = [5:7:7*Neb-2].';
   i7f = [6:7:7*Neb-1].';
   i7g = [7:7:7*Neb].';

   bst = [i7a;i7b;i7c;i7d;i7e;i7f;i7g];

   i7a = [1:7:7*Ncp-6].';
   i7b = [2:7:7*Ncp-5].';
   i7c = [3:7:7*Ncp-4].';
   i7d = [4:7:7*Ncp-3].';
   i7e = [5:7:7*Ncp-2].';
   i7f = [6:7:7*Ncp-1].';
   i7g = [7:7:7*Ncp].';

   sind = [i7a;i7b;i7c;i7d;i7e;i7f;i7g];

   if (a.icp(1) < 0)

      % Modal aero DOF reduction.
      Psix = bladeModeAero (bst,sind,bsh,-a.icp);

   elseif (a.icp(1) > 0)

      % Spline aero DOF reduction.
      Psix = splineAero (bst,sind,a.icp,rp(1:Neb));

   end

   Neb7 = 7*Neb;
   Ncp7 = 7*Ncp;
   Psi  = [Psix sparse(Neb7,Ncp7) sparse(Neb7,Ncp7); ...
           sparse(Neb7,Ncp7) Psix sparse(Neb7,Ncp7); ...
           sparse(Neb7,Ncp7) sparse(Neb7,Ncp7) Psix];

end