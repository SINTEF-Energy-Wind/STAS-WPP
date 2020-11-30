function [M,G,H,C,K,Q,dM,dMg,dG,dGd,dH,dHd,dC,dCd,dK,dQ,dQd] = ...
                     buildMRQ (linFlag,s,q,qd,qdd,P,F,grav)
%
% Assemble the mass matrix and force vectors, and their linearizations.
% This basically converts G and C from matrices to force vectors,
% by multiplying with qd, and adding soil reactions and other applied
% force effects via the Qh matrix.
%
% Version:        Changes:
% --------        -------------
% 16.01.2020      Original code, based on buildMRQLin/NL.
%
% Version:        Verification:
% --------        -------------
% 16.01.2020      Derivatives verified by complex step on a selection
%                 of reference, joint, and elastic DOFs.
%
% Inputs:
% -------
% linFlag         : = 1 for linearization.
% s               : Wind turbine data structure.
% q,qd,qdd        : Nodal positions, velocities, accelerations in the
%                   "q" basis.
% P               : Nodal positions.
% F               : Initial forces.
% shape           : Body modes.
% grav            : Gravitational vector in global coordinates.
%
% Outputs:
% --------
% M               : Mass matrix.
% G               : Vector, Gmat-times-dq/dt.
% H               : H vector.
% C               : Vector, Cmat-times-dq/dt.
% K               : K vector.
% dM              : dM/dq*d2q/dt2.
% dMg             : dM/dq*gravitational acceleration.
% dG              : dG/dq
% dGd             : dG/dqd
% dH              : dH/dq
% dHd             : dH/dqd
% dC              : dC/dq
% dCd             : dC/dqd
% dK              : dK/dq
% T               : Transform of basis [x;xd] = (T_z^x)*[z;zd].
% dT              : (dT/dq)*[qd;qdd]
% dTd             : (dT/dqd)*[qd;qdd]

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
Ndj  = size(P,1);

% Call to construct the nonlinear terms.  Recall that only the
% nonlinear outputs are valid.
[Mf,dMf,dMgf,Gf,dGf,dGfd,Hf,dHf,dHfd,Cf,dCf,Kf,dKf] = ...
                       assembleElements (0,s,q,qd,qdd,P,grav);

% Add soil forces.
Ndf = 6*s.foundation.Nnod;
dofs = idofs(1) + [1:Ndf].';
Fsoil = soilNL (s.foundation.k(1,:).',s.foundation.k(2,:).',               ...
        s.foundation.k(3,:).',s.foundation.k(4,:).',s.foundation.k(5,:).', ...
        s.foundation.k(6,:).',s.foundation.k(7,:).',s.foundation.k(8,:).', ...
        P(dofs),q(dofs),qd(dofs));
Fnod = F;
Fnod(dofs) = Fnod(dofs) + Fsoil;

% Q matrix for applied forces.
Qh = Qhat (q,P,idofs,inods);

% RHS forces, unconstrained.  Nonlinear outputs.  R = -G + H - C - K + Q.
M = Mf;
G = Gf*qd;
H = Hf;
C = Cf*qd;
K = Kf;
Q = Qh*Fnod;

if (linFlag == 1)

   % Assemble the forces and linearizations.  Everything is in place,
   % like dG = dG/dq*dq/dt, except the G vector and C vector.  For
   % these the matrix needs to be multiplied with the dq/dt vector
   % and the derivatives wrt qd augmented.
   [M,dM,dMg,Gf,dG,dGfd,H,dH,dHd,Cf,dC,K,dK] = ...
                       assembleElements (1,s,q,qd,qdd,P,grav);

   G = Gf*qd;
%  dG = dGf;            % (dG/dq)*qd, direct from assembleElements.
   dGd = Gf + dGfd;     % d/dqd (Gf*qd) = Gf + (dGf/dqd)*qd.

   C = Cf*qd;
   dCd = Cf;

   % Now the applied forces.  Prepare soil forces and derivatives.
   [Fsoil,dFsoildq,dFsoildqd] =                                            ...
        soilLin (s.foundation.k(1,:).',s.foundation.k(2,:).',              ...
        s.foundation.k(3,:).',s.foundation.k(4,:).',s.foundation.k(5,:).', ...
        s.foundation.k(6,:).',s.foundation.k(7,:).',s.foundation.k(8,:).', ...
        P(dofs),q(dofs),qd(dofs));

   % Mean forces acting through displacements.
   dQ = dQhdqF (q,P,Fnod,idofs,inods); 

   % Fluctuating soil forces.
   dQ(:,dofs) = dQ(:,dofs) + Qh(:,dofs)*dFsoildq;

   dQd = sparse (Ndj,Ndj);
   dQd(:,dofs) = Qh(:,dofs)*dFsoildqd;

else

   dM  = sparse(Ndj,Ndj);
   dMg = sparse(Ndj,Ndj);
   dG  = sparse(Ndj,Ndj);
   dGd = sparse(Ndj,Ndj);
   dH  = sparse(Ndj,Ndj);
   dHd = sparse(Ndj,Ndj);
   dC  = sparse(Ndj,Ndj);
   dCd = sparse(Ndj,Ndj);
   dK  = sparse(Ndj,Ndj);
   dQ  = sparse(Ndj,Ndj);
   dQd = sparse(Ndj,Ndj);

end
