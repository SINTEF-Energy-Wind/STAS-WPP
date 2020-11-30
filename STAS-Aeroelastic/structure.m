function [Lmat,Rvec,y,A,Bu,By,C,Du,Dy,ret,slv] =  ...
               structure (linFlag,Tflag,s,x,dx,P,F,shape,mdamp,grav)
%
% Assemble the structural equations in state space form.  Note that
% q = (Lamq shape) eta, where Lamq has undone the partitioning of 
% Lambda.
%
%   States:              y vector:             u vector:
%   eta      Neta        q          Ndj        F          Ndj
%   deta/dt  Neta        dqdt       Ndj        d2eta/dt2  Neta
%                        d2qdt2     Ndj
%
% Version:        Changes:
% --------        -------------
% 22.05.2018      Original code.
% 18.01.2020      Rewritten to accommodate updated buildMRQ and
%                 transforms of the basis.
%
% Version:        Verification:
% --------        -------------
% 22.05.2018      Match with structureNL verified by complex step.
%                 Converged Newton solutions for nonlinear equations.
%                 NREL 5 MW turbine modes match published values.
% 18.01.2020      A-matrix derivatives of Lmat and Rvec verified by
%                 complex step.  C-matrix derivatives of y verified
%                 by complex step.
%
% Inputs:
% -------
% linFlag         : = 1 for linearization.
% Tflag           : = 0 if the basis is [q,qd]
%                   = 1 if the basis is [qh,qhd]
%                   = 2 if the basis is [eta,etad]
% s               : Data structure describing the structure.
% x               : [eta;etad].
% dx              : d/dt [eta;etad].
% P               : Undeformed nodal positions.
% F               : Ndj nodal forces.
% shape,mdamp     : Mode shapes and modal damping vector.
% grav            : Gravity vector in global coordinates.
%
% Outputs:
% --------
% Lmat,A...Dy     : State matrices.
% Rvec            : NL RHS vector of forces.
% ret,slv         : Retained and slave DOFs.

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

Nx   = size(x,1);
Neta = Nx/2;
Ndj  = size(P,1);
Nu   = Ndj + Neta;

iep = [1:Neta].';
iev  = Neta + [1:Neta].';

iqp = [1:Ndj].';
iqv = Ndj + [1:Ndj].';

% Get retained and slave DOFs.
slv = slaveDOFs (idofs);
vec = [1:Ndj].';
[jnk,ret,jnk2] = partitionMatrix (vec,slv,[]);

% In order to call buildMRQ, which is in the q basis, we need to
% get there from the input.  This is nominally the eta basis
% (constrained, and reduced with body modes) but this can be
% bypassed using the Tflag input.
%
% Build the nonlinear transforms, without knowing the acceleration.
eta   = x(iep);
etad  = x(iev);
etadd = dx(iev);

if (Tflag >= 2)
   [qh,qhd,qhdd,Tetaqh] = buildTetaqh (eta,etad,etadd,shape);
else
   qh     = eta;
   qhd    = etad;
   qhdd   = etadd;
   Tetaqh = speye(2*size(qh,1));
end

if (Tflag >= 1)
   [q,qd,qdd,Tqhq,jnk,jnk2] = buildTqhq (0,s,qh,qhd,qhdd,P,ret,slv);
else
   q    = qh;
   qd   = qhd;
   qdd  = qhdd;
   Tqhq = speye(2*size(q,1));
end

% Overall transform matrix.
T = Tqhq*Tetaqh;
TT = T.';

% Build the nonlinear terms.  Tflag = 2, basis is eta.
[M,G,H,C,K,Q,dM,dMg,dG,dGd,dH,dHd,dC,dCd,dK,dQ,dQd] = ...
                      buildMRQ (0,s,q,qd,qdd,P,F,grav);

y = [q;qd;qdd];

% This should reproduce the effective gravity vector implemented
% by assembleElements when computing dMg.  Each body reference
% node is given an acceleration of g.
gvec = zeros(Ndj,1);
indb = [1 2 3 4 6 7 8].';  % Joint 5 is not a body reference.
for ibod = 1:7
   gvec(idofs(indb(ibod))+[1:3]) = grav;
end
Mg = M*gvec;

% Build and transform/reduce the state equations for the structure.
Rq = -G+H-C-K+Q+Mg;

Lmat = sparse(Nx,Nx);
Lmat(iep,iep) = speye(Neta);
Lmat(iev,:) = real(TT(iev,iqv))*M*T(iqv,:);  % Shp'*Lam'*M*[Gam*Shp Lam*Shp]

Rvec = zeros(Nx,1);
Rvec(iep) = x(iev);
Rvec(iev) = real(TT(iev,iqv))*Rq;

if (real(s.zdamp) > 0)
   Rvec(Neta+[1:Neta]) = Rvec(Neta+[1:Neta]) ...
                       - mdamp.*x(Neta+[1:Neta]);
end

if (linFlag == 1)

   etadd = dx(Neta+[1:Neta]);

   % Get the linearized transforms.
   if (Tflag >= 2)
      [qh,qhd,qhdd,Tetaqh] = buildTetaqh (eta,etad,etadd,shape);
   else
      qh     = eta;
      qhd    = etad;
      qhdd   = etadd;
      Tetaqh = speye(2*size(qh,1));
   end

   if (Tflag >= 1)
      [q,qd,qdd,Tqhq,dT,dTd] = buildTqhq (1,s,qh,qhd,qhdd,P,ret,slv);
   else
      q    = qh;
      qd   = qhd;
      qdd  = qhdd;
      Tqhq = speye(2*Ndj);
      dT   = sparse(2*Ndj,Ndj);
      dTd  = sparse(2*Ndj,Ndj);
   end

   % (Note, dQd is not due to the Qhat matrix, which is only a function
   %  of q, but rather it is due to the soil damping.)
   [M,G,H,C,K,Q,dM,dMg,dG,dGd,dH,dHd,dC,dCd,dK,dQ,dQd] = ...
                     buildMRQ (1,s,q,qd,qdd,P,F,grav);
   Qh = Qhat (q,P,idofs,inods);

   dRq  = [-dG+dH-dC-dK+dQ+dMg, -dGd+dHd-dCd+dQd];
   dRdamp = [sparse(Neta,Neta), diag(mdamp)];
   AM = dM;
   AT = M*dT(iqv,:);
   ATd = M*dTd(iqv,:);

   A = sparse(Nx,Nx);
   A(iep,iev) = speye(Neta);
   A(iev,:) = real(TT(iev,iqv))*(dRq - [AM+AT, ATd])*T - dRdamp;
   Bu = sparse(Nx,Nu);
   Bu(iev,:) = real(TT(iev,iqv))*[Qh, sparse(Ndj,Neta)];
   By = sparse(Nx,3*Ndj);
   C  = [T;                               ...
         [sparse(Ndj,Neta), T(iqv,1:Neta)] ...
      +   dT(iqv,:)*T(1:Ndj,:)             ...
      +   dTd(iqv,:)*T(iqv,:)];
   Du = [sparse(2*Ndj,Ndj+Neta); sparse(Ndj,Ndj), T(iqv,Neta+[1:Neta])];
   Dy = sparse(3*Ndj,3*Ndj);

else

   A = sparse(Nx,Nx);
   Bu = sparse(Nx,Ndj+Neta);
   By = sparse(Nx,3*Ndj);
   C  = sparse(3*Ndj,Nx);
   Du = sparse(3*Ndj,Ndj+Neta);
   Dy = sparse(3*Ndj,3*Ndj);

end