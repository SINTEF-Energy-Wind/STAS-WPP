function Psix = bladeModeAero (bst,mind,bshape,imd)
% 
% This function defines a basis for the reduction of aerodynamic
% states from the selected structural mode shapes.
%
% Version:        Changes:
% --------        -------------
% 16.03.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 16.03.2018      Checked with sample inputs.
%
% Inputs:
% -------
% bst             : A vector of indices, Neb*Nxe in size, containing
%                   the positions of each group of Neb aero states
%                   that are to be reduced according to the basis
%                   functions.
% mind            : A vector of indices, Nimd*Nxe in size, 
%                   specifying the order of modal states.
% bshape          : The body mode shapes associated with the blade
%                   elements.  Neb rows by Nmod columns.
% imd             : The modes, out of the Nmod input in bshape, to use
%                   as basis functions.
%
% Outputs:
% --------
% Psix            : The matrix used to reduce the aero states.

Neb  = size(bshape,1);
Nx   = size(bst,1);
Nxe  = Nx/Neb;
Nimd = size(imd,1);
Nz   = Nimd*Nxe;

Psix = spalloc (Nx,Nz,Nx*Nimd);
for ixe = 1:Nxe
   icb =  Neb*(ixe-1);
   ici = Nimd*(ixe-1);
   Psix(bst(icb+[1:Neb]),mind(ici+[1:Nimd])) = bshape(:,imd);
end

%[Neb Nx Nxe Nimd Nz]
