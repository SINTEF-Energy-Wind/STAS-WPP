function [dxadt,ya] = aeroNL (psiFlag,s,a,xa,t,q,dqdt,P, ...
                              Tas,Try,Vg,ch,Lel,foilwt,  ...
                              aoaz,aoast,xas,yas,Psi)
%
% Prepares and executes the call to BEMNL for each annulus, computing
% quantities that are updated on each timestep.
%
%   States:           y vector:
%   ad        1       Fa        1:6    (Repeat for each blade element.)
%   a1,a2    2:3      
%   Vih z,t  4:5      
%   Vi z,t   6:7
%
% Version:        Changes:
% --------        -------------
% 16.04.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 16.04.2018      
%
% Inputs:
% -------
% psiFlag         : Set to 1 to implement the dynamic wake in multi-blade
%                   coordinates.  Set to zero for blade-by-blade. 
%                   Regardless, dxdt is reported as blade-by-blade.
%                   If psiFlag = 1, then it is assumed that the elements
%                   are packed as Bl1 e1 e2 e3 ... eN, Bl2 e1 e2 ... for
%                   three blades.
% s,a             : Structural and aerodynamic input data structures.
% xa              : The vector of aero states.  For each element, these
%                   are ad, a1, a2, Virhz, Virht, Virz, and Virt. 
% t               : Time.
% q,dqdt          : Full structural displacement and velocity DOFs.
% P               : Nodal positions.
% Tas             : Transform from airfoil to section coordinates.
% Try             : Transform from rotorplane to yaw coordinates.
% Vg              : 3*Nel vector, incoming windspeed x,y,z in global
%                   coordinates.
% ch              : Airfoil chord length.
% Lel             : Blade element length.
% foilwt          : Nfoil-by-Nel table.  Weights to use when computing
%                   airfoil coefficients from splined tables.
% aoaz            : Zero-lift angles-of-attack for each element. Should
%                   be computed precisely from the airfoil tables.
% aoast           : 2*Nel vector, containing deep-stall angles-of-attack
%                   for each element.  Alternating positive, negative.
% xas,yas         : X^a and Y^a coordinates of the reference section
%                   coordinate system.
% Psi             : Matrix of basis functions for aero states.
%
% Outputs:
% --------
% dxadt           : Rate of change of states xa.
% ya              : Vector of outputs which are to be passed to other
%                   modules.

%'aeroNL'

Nel = a.Nb*a.Neb;
Neb = a.Neb;
Ndj = size(q,1);
Nxf  = 7*Nel;

xaf = Psi*xa;

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
Ydof  = idofs(3);
Ddof  = idofs(4);
nodof = idofm(6) - 6; % idofs(5);
[omega,azi,jnk] = rotorSpeedAero (q,dqdt,P,Try,Ydof,Ddof,nodof);

[Tar,Trg,TB0g,TsB,TBB0,dTar,dTsB,dTBB0,wg] = ...
                 BEMprepTransforms (s,a,q,dqdt,P,Tas);
[zr,Area,Dp,r,Lp,xeg,xhg,xyg] = ...
                 BEMprepProjections (s,a,q,P,Try,Trg);

dxafdt = zeros(Nxf,1);
ya    = zeros(6*Nel,1);
for iel = 1:Neb  % Could use parallel processing here.

%iel

   ind = [iel Neb+iel 2*Neb+iel].';

   i7 = 7*(ind-1);
   i6 = 6*(ind-1);
   i3 = 3*(ind-1);
   i2 = 2*(ind-1);

   ind7 = [i7(1)+[1:7] i7(2)+[1:7] i7(3)+[1:7]].';
   ind6 = [i6(1)+[1:6] i6(2)+[1:6] i6(3)+[1:6]].';
   ind3 = [i3(1)+[1:3] i3(2)+[1:3] i3(3)+[1:3]].';
   ind2 = [i2(1)+[1:2] i2(2)+[1:2] i2(3)+[1:2]].';

   [dxafdt(ind7),ya(ind6)] = BEMNL (psiFlag,                ...
                             xaf(ind7),t,Tar(:,ind3),Trg,   ...
                             Vg(ind3),wg(ind3),             ...
                             zr(ind2),ch(ind),Lel(ind),     ...
                             a.aoas,a.kfoils,foilwt(:,ind), ...
                             aoaz(ind),aoast(ind2),         ...
                             xas(ind),yas(ind),a.dens,      ...
                             Area(ind),Dp,azi,omega);

end


dxadt = (Psi.')*dxafdt;
