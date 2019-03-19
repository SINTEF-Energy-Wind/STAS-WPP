function [Lmat,Rvec,y] =                                ...
               aeroelasticNL (t,x1,y1,Fext,grav,        ...
                              psiFlag,s,a,P,            ...
                              ret,slv,shape,mdamp,      ...
                              Tas,Try,Vg,ch,Lel,foilwt, ...
                              aoaz,aoast,xas,yas,Psi)
%
% Build the nonlinear aeroelastic state equations at a specified
% point in time, where a state and output vector are available as
% input.
%
%   States:           y vector:
%   eta       N       q            Ndj
%   deta/dt   N       dqdt         Ndj
%                     d2qdt2       Ndj
% -----
% (For each aero element)
%   ad         1      Fa            6
%   a1,a2     2:3
%   Vih z,t   4:5
%   Vi z,t    6:7
% -----
%                     F            Ndj
%
% Version:        Changes:
% --------        -------------
% 17.04.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 17.04.2018      
%
% Inputs:
% -------
% t               : Time (presently unused).
% x1,y1           : Present states and associated outputs (only q
%                   and dqdt are used).
% Fext            : Additional nodal forces that are not part of
%                   aeroNL.
% grav            : Gravitational acceleration in global CS.
% psiFlag         : = 1 for MBC dynamic wake, 0 for independent-
%                   blade dynamic wake.
% s,a             : Structural and aero data structures.
% P               : Nodal offsets.
% ret,slv         : Lists of retained and slave DOFs.
% shape           : Mode shape matrix.
% mdamp           : Modal damping vector.
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
% 

%'---------------------'
%'NL'

Nxs = 2*size(shape,2);
Nxa = size(Psi,2);
Nx  = size(x1,1);
Ny  = size(y1,1);
Ndj = size(P,1);

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

xs = x1(1:Nxs);
xa = x1(Nxs+[1:Nxa]);

q    = y1(1:Ndj);
dqdt = y1(Ndj+[1:Ndj]);

% Aerodynamic state equations.
[dxadt,ya] = aeroNL (psiFlag,s,a,xa,t,q,dqdt,P, ...
                     Tas,Try,Vg,ch,Lel,foilwt,  ...
                     aoaz,aoast,xas,yas,Psi);

% Update force vector.
Fa = ya;
F = airfoilToNodal (s,a,Fa,Tas,q,P,idofs);
F = F + Fext;

% Structural state equations.
[Lms,Rvs,ys] = structureNL (xs,t,s,P,F,shape,mdamp,grav, ...
                            ret,slv,q(slv));

Lmat = [Lms sparse(Nxs,Nxa);sparse(Nxa,Nxs) speye(Nxa)];
Rvec = [Rvs;dxadt];
y    = [ys;ya;F];


