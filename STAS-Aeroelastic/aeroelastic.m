function [Lmat,Rvec,y,A,Bu,By,C,Du,Dy] =           ...
            aeroelastic (linFlag,s,a,x,u,P,        ...
                         shape,mdamp,grav,         ...
                         Tas,Try,ch,Lel,foilwt,    ...
                         aoaz,aoast,xas,yas,Psi)
%
% Compute the terms in the nonlinear aeroelastic equations and build a
% linear model of the turbine at a given operating point, which does
% not need to be an equilibrium point.  In the output matrices, joints
% are free, and need to be controlled, or locked by eliminating their
% DOFs.
%
% The interface (input and output) variables for the nonlinear model
% include nodal forces F on the structure; wind speeds Vg, in global
% coordinates, at each blade element; and nodal displacements,
% velocities, and accelerations in the "q" basis.  The linear model
% includes many additional quantities that are interface variables
% between the structural and aerodynamic models, as well as between
% modules within the aerodynamic model.
%
% ================
% Nonlinear model:
%   States:              y vector:             u vector:
%   eta       Neta       q            Ndj      F          Ndj
%   deta/dt   Neta       dqdt         Ndj      d2eta/dt2  Neta
%                        d2qdt2       Ndj
% -----
% (For each aero element)
%   ad         1         Fa            6       Vg          3
%   a1,a2     2:3
%   Vih z,t   4:5
%   Vi z,t    6:7
% -----
%                        F            Ndj
%
% ================
% Linear model, intermediate variables required for automated linking
% in linkAeroelastic.
%   States:              y vector:             u vector:
% ----------------------- Structure --------------------------
%   eta       Neta       q         Ndj         F          Ndj
%   deta/dt   Neta       dq/dt     Ndj         d2eta/dt2  Neta
%                        d2q/dt2   Ndj
%                        xng     3*Nnod
%                        vng     6*Nnod
%                        F         Ndj
% ---------------------- Aerodynamic -------------------------
%   ad         1         (q)       1:24        Vg          3
%   a1,a2     2:3        (dq/dt)  25:48
%   Vih z,t   4:5        xng1,2   49:54
%   Vi z,t    6:7        vng1,2   55:60
%                        wg       61:63
%                        wa       64:66
%                        aq        67
%                        Cl,Cd,Cm 68:70
%                        Fl,Fd,M  71:73
%                        Fa       74:79
%                        Fp       80:85
%                        Fr       86:91
%                        Fzts     92:94
%                        Vg       95:97
%                        Ua       98:100
%                        Umag      101
%                        Ur      102:104
%                        Uzts    105:107
%                        Vzts    108:110
%                        Viq     111:112
%                        Viy     113:114
%                        Vixyz   115:117
%                        Wmag      118
%                        xeg     119:121
%                        xhg     122:124
%                        xnr1,2  125:130  
%                        xer     131:133
%                        r         134
%                        Lp        135
%                        z         136
%                        f         137  
% (Repeat the above x,y,u consecutively for each blade element.)
%                        Azi       (1)  (effective rotor azimuth)
%                        Waero     (1)  (aero rotor speed)
%                        Dp        (3)
%
%
% Version:        Changes:
% --------        -------------
% 06.02.2018      Original code.
% 21.01.2020      Revised to accommodate the updates to structure.m.
%
% Version:        Verification:
% --------        -------------
% 06.02.2018      I expect there is a problem with the (dM/dq)(d^2q/dt^2)
%                 term's linearization, but it is tricky to find, and it
%                 is generally not critical, as the linearization is
%                 seldom relevant in a condition with high steady DOF
%                 accelerations.
% 21.01.2020      "A" matrix derivatives show precise agreement, based on
%                 complex step using Lmat, Rvec.  The above issue with
%                 accelerations is no longer evident.
%
% Inputs:
% -------
% linFlag         : = 1 for linearization.
% s               : Data structure describing the structures.
% a               : Data structure describing the aerodynamics.
% x               : State vector.
% u               : Input vector.
% P               : Undeformed nodal positions.
% shape,mdamp     : Mode shapes and modal damping vector.
% grav            : Gravity vector in global coordinates.
% Psi             : Aero mode shapes.
%
% Outputs:
% --------
% Lmat, Rvec      : LHS and RHS of nonlinear state equation.
% y               : Nonlinear y vector.
% A...Dy          : State matrices.

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

Neta = size(shape,2);
Nxs  = 2*Neta;
Nxa  = size(Psi,2);
Ndj  = size(P,1);
Nel  = s.blade(1).Nel + s.blade(2).Nel + s.blade(3).Nel;
Nnod = Ndof/6;

Nx   = size(x,1);
Nu   = Ndj + Neta + 3*Nel;
Ny   = 4*Ndj + 137*Nel + 9*Nnod + 5;

% Retained and slave structural DOFs.
slv = slaveDOFs (idofs);
vec = [1:Ndj].';
[jnk,ret,jnk2] = partitionMatrix (vec,slv,[]);

% Compute q DOFs from x input.
xs   = x(1:Nxs);
eta  = xs(1:Neta);
etad = xs(Neta+[1:Neta]);
zro  = zeros(Neta,1);     % Not needed for nonlinear.
[qh,qhd,jnk1,jnk2] = buildTetaqh (eta,etad,zro,shape);
[q,qd,jnk2,jnk3,jnk4,jnk5] = buildTqhq (0,s,qh,qhd,jnk1,P,ret,slv);

% Aerodynamic state equations.  Hard-code psiFlag = 1, wake dynamics
% in MBC coordinates.
xa = x(Nxs+[1:Nxa]);
Vg = u(Ndj+Neta+[1:3*Nel]);
[dxadt,ya] = aeroNL (1,s,a,xa,0,q,qd,P,        ...
                     Tas,Try,Vg,ch,Lel,foilwt, ...
                     aoaz,aoast,xas,yas,Psi);

% Update force vector.
Fext = u(1:Ndj);
F = airfoilToNodal (s,a,ya,Tas,q,P,idofs);  % Fa = ya.
F = F + Fext;

% Structural state equations.
dxs  = zeros(Nxs,1);
dxs(1:Neta) = xs(Neta+[1:Neta]);
dxs(Neta+[1:Neta]) = u(Ndj+[1:Neta]);
[Lms,Rvs,ys,jA,jBu,jBy,jC,jDu,jDy,ret,slv] =  ...
               structure (0,2,s,xs,dxs,P,F,shape,mdamp,grav);

Lmat = [Lms sparse(Nxs,Nxa);sparse(Nxa,Nxs) speye(Nxa)];
Rvec = [Rvs;dxadt];
y    = [ys;ya;F];   % F output: total forces.

if (linFlag == 1)

   % =========== Structural state equations =============
   [lls,Rvs,ys,aas,bbus,bbys,ccs,ddus,ddys,ret,slv] =  ...
             structure (1,2,s,xs,dxs,P,F,shape,mdamp,grav);
   Nxs  = size(aas,1);
   Nus  = size(bbus,2);
   Nysq = size(ccs,1);
   Nys  = Nysq + 9*Nnod + Ndj;
   Nxs2 = Nxs/2;

   % =========== Aerodynamic state equations ============ 
   if (a.icp(1) == 0)
      Ncp = a.Neb;
   else
      Ncp = size(a.icp,1);
   end
   bsh = bladeModeShape (s,ret,shape);
   [dxadt,aaa,bbya,cca,ddya,jnk] = ...
              aeroLin (1,xa,Vg,s,a,q,qd,P,bsh);
   Nxa = size(aaa,1);
   Nya = size(cca,1);
   Nua = 3*a.Nb*a.Neb;  % We will define Vg as a global input.
   Nyae = 137;          % Number of aero y's associated with each element.
   Nea = a.Nb*a.Neb;
   Neb = a.Neb;
   lla = speye (Nxa);

   % ======= Equations for linking aero-structure ======= 
   [Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);
   [xng,vng,ddypv] = globalPosVel (s,idofs,q(1:Ndof),P(1:Ndof),qd(1:Ndof));
   [Waero,azi,ddyRSA] = rotorSpeedAero (q,qd,P,Tn_y, ...
                                        idofs(3),idofs(4),idofm(6)-6);

   [llas,aaas,bbuas,bbyas,ccas,dduas,ddyas] =                       ...
                  linkAeroelastic (s,a,                             ...
                                   lls,aas,bbus,bbys,ccs,ddus,ddys, ...
                                   lla,aaa,bbya,cca,ddya,           ...
                                   ddypv,ddyRSA);

   A    = aaas;
   Bu   = bbuas;
   By   = sparse(bbyas);
   C    = sparse(ccas);
   Du   = sparse(dduas);
   Dy   = sparse(ddyas);

else

   A  = sparse(Nx,Nx);
   Bu = sparse(Nx,Nu);
   By = sparse(Nx,Ny);
   C  = sparse(Ny,Nx);
   Du = sparse(Ny,Nu);
   Dy = sparse(Ny,Ny);

end




