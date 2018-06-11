function Psix = splineAero (bst,sind,icp,rel)
%
% Splines are used to reduce the number of aerodynamic states.
% This function produces the matrix used to perform the
% reduction of the aero states.
%
% Version:        Changes:
% --------        -------------
% 16.03.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 16.03.2018      
%
% Inputs:
% -------
% bst             : A vector of indices, Neb*Nxe in size, containing
%                   the positions of each group of Neb aero states
%                   that are to be reduced according to the basis
%                   functions.
% sind            : A vector of indices, Ncp*Nxe in size, specifying 
%                   the order of spline states.
% icp             : A list of elements to use as control points
%                   for the spline curves.
% rel             : The radial coordinate of each blade element.
%
% Outputs:
% --------
% Psix            : The matrix used to reduce the aero states.

Neb = size(rel,1);

rcp = rel(icp);
Ncp = size(icp,1);

Nx  = size(bst,1);
Nxe = Nx/Neb;
Nz  = Ncp*Nxe;

[F,G] = splineCoefficients (rcp);
Sp = splineSMatrix (rcp,rel);

FiG = F\G;
H = Sp*FiG;  % The elemental transformation for each group of
             % Neb states.

Psix = spalloc (Nx,Nz,Nx*Ncp);
for ixe = 1:Nxe
   icb = Neb*(ixe-1);
   icc = Ncp*(ixe-1);
   Psix(bst(icb+[1:Neb]),sind(icc+[1:Ncp])) = H;
end
