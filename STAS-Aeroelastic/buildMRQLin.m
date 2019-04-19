function [M,dM,MG,dMG,dMGd,R,dR,dRd,Q,dQ,dQd,slv,ret, ...
          Lambda,Gamma,Leq,dLdq,dLam,dGam,dGamd] =    ...
                         buildMRQLin (s,q,dqdt,d2qdt2,P,F)
%
% Assemble the linearized mass, gyroscopic/damping, and stiffness
% matrices for the turbine structure.  Note that
% R = G-H+C+K-L+Q from buildMGHKfull.
%
% Version:        Changes:
% --------        -------------
% 23.11.2017      Original code.
% 18.05.2018      Updated with assembleElementsLin as a direct
%                 linearization of assembleElementsNL.
% 24.02.2019      Deleted call to hydrodynamic added mass function;
%                 now this is handled as extra mass in the element
%                 mass matrices.
%
% Version:        Verification:
% --------        -------------
% 23.11.2017      M,C,K checked against buildMGHKfull, using complex
%                 step.
% 18.05.2018      
% 24.02.2019      
%
% Inputs:
% -------
% s               : Wind turbine data structure.
% q...d2qdt2      : Initial displacements, velocities, accelerations.
% P               : Nodal positions.
% F               : Initial forces.
%
% Outputs:
% --------
% M               : Mass matrix.
% dM              : Matrix dM/dq*d2q/dt2.
% MG ... dMGd     : M*Gamma*dqh/dt terms.
% R               : vector, = G-H+C+K-L, constrained via Lambda^T.
% dR              : dR/dqh matrix, derivatives wrt reduced qhat DOFs.
%                   The negative of the stiffness, in other words.
% dRd             : dR/dqhdot matrix.  Damping, centrifugal, and
%                   gyroscopic terms.  The damping includes Rayleigh
%                   damping, if this is active, but not modal damping.
% Q               : Vector of generalized forces.
% dQ              : dQ/dqh.
% dQd             : dQ/dqhdot.
% slv,ret         : Lists of slave and retained DOFs.
% Lambda          : DOF constraint matrix.
% Gamma           : Constraint force matrix.
% Leq             : The L matrix, constraint equation Jacobian.
% dLdq            : dL/dq, second partial derivatives of the constraint
%                   equations.
% dL              : dLambda/dq*d2q/dt2.
% dGam            : dGamma/dq*dq/dt.
% dGamd           : dGamma/dqdot*dq/dt.

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);

Ndj = size(q,1);

% Compute constraint matrices Lambda, containing the partitioned
% constraint equations, and Gamma, a term involved in computing the
% constraint forces.
gdofs = getgdofs (idofs,idofm,Ndj);
[Lambda,Leq,Con,ret,slv] = constraints (q,P,Tb_h,idofs,idofm);
[Gamma,dLdq] = constraintGamma (q,dqdt,P,Leq,Lambda, ...
                                Tb_h,ret,slv,idofs,idofm);
Nret = size(ret,1);
Nslv = size(slv,1);

LamT = Lambda.';

% Build the linearized terms, at this point unconstrained and
% with Rayleigh damping.  (Modal damping, if active, comes later.)
rtsl = [ret;slv];
d2qL = zeros(Ndj,1);
d2qG = zeros(Ndj,1);
d2qL(rtsl) = Lambda*d2qdt2(ret);
d2qG(rtsl) = Gamma*dqdt(ret);
[M1,dML,dMG,G1,dG1,dGd1,H1,dH1,dHd1,C1,dC1,K1,dK1] =       ...
               assembleElementsLin (s,idofs,q,dqdt,d2qL,d2qG,P);

[M1p,rr,cr]   = partitionMatrix (M1,slv,slv);
clear M1;
[dMLp,rr,cr]  = partitionMatrix (dML,slv,slv);
clear dML;
[dMGp,rr,cr]  = partitionMatrix (dMG,slv,slv);
clear dMG;
[G1p,rr,cr]   = partitionMatrix (G1,slv,slv);
clear G1;
[dG1p,rr,cr]  = partitionMatrix (dG1,slv,slv);
clear dG1;
[dGd1p,rr,cr] = partitionMatrix (dGd1,slv,slv);
clear dGd1;
[H1p,rr,cr]   = partitionMatrix (H1,slv,[]);
clear H1;
[dH1p,rr,cr]  = partitionMatrix (dH1,slv,slv);
clear dH1;
[dHd1p,rr,cr] = partitionMatrix (dHd1,slv,slv);
clear dHd1;
[C1p,rr,cr]   = partitionMatrix (C1,slv,slv);
clear C1;
[dC1p,rr,cr]  = partitionMatrix (dC1,slv,slv);
clear dC1;
[K1p,rr,cr]   = partitionMatrix (K1,slv,[]);
clear K1;
[dK1p,rr,cr]  = partitionMatrix (dK1,slv,slv);
clear dK1;

[qp,rr,cr]    = partitionMatrix (q,slv,[]);
[dqdtp,rr,cr] = partitionMatrix (dqdt,slv,[]);
[d2qdt2p,rr,cr] = partitionMatrix (d2qdt2,slv,[]);

% ----------------------------------------------------------------
% Add soil stiffness and damping to the foundation.  Consider the
% stiffness.  Each nodal spring applies a force F = -ks*d, with d
% the global lateral displacement.  Together, these forces can be
% written F = -Ks*q.  Applying the nodal forces requires the
% operation Qh*F.  Moving to the LHS, the effective stiffness is
% Qh*Ks.
Qh = Qhat (q,P,idofs,inods);
Ndf = 6*s.foundation.Nnod;
dofs = idofs(1) + [1:Ndf].';
Qhf = Qh(dofs,dofs);

[Fsoil,dFsoildq,dFsoildqd] =                                               ...
        soilLin (s.foundation.k(1,:).',s.foundation.k(2,:).',              ...
        s.foundation.k(3,:).',s.foundation.k(4,:).',s.foundation.k(5,:).', ...
        s.foundation.k(6,:).',s.foundation.k(7,:).',s.foundation.k(8,:).', ...
        P(dofs),q(dofs),dqdt(dofs));

% Add the effect of soil forces.
Fnod = F;
Fnod(dofs) = Fnod(dofs) + Fsoil;

% ----------------------------------------------------------------
% Generalized forces.
QF = Qh*Fnod;
[QFp,rr,cr] = partitionMatrix (QF,slv,[]);
Q = LamT*QFp;



%full([[1:Ndj].' Fnod QF])
%full([[1:Nret].' Q])



% Mean forces acting through displacements.
dQfull = dQhdqF (q,P,Fnod,idofs,inods); 

% Force fluctuations.  
dQfull(:,dofs) = dQfull(:,dofs) + Qh(:,dofs)*dFsoildq;

[dQp,rr,cr] = partitionMatrix (dQfull,slv,slv);



%matp = dQp*Lambda;
%rtsl = [ret;slv];
%mat = sparse(size(matp,1),size(matp,2));
%mat(rtsl,:) = matp;
%full([[1:Ndj].' mat(:,78)])




dQdfull = sparse (Ndj,Ndj);
dQdfull(:,dofs) = Qh(:,dofs)*dFsoildqd;
[dQdp,rr,cr] = partitionMatrix (dQdfull,slv,slv);

dLTdqQ = dLambdaTdq (Lambda,Leq,dLdq,gdofs,ret,slv,QFp);
dQ = LamT*(dQp*Lambda + dQdp*Gamma) + dLTdqQ*Lambda;

dQd = LamT*dQdp*Lambda;

% ----------------------------------------------------------------
% Mass matrix.
M  = LamT*M1p*Lambda;
MG = LamT*M1p*Gamma;

dLTdqMLd2q = dLambdaTdq (Lambda,Leq,dLdq,gdofs,ret,slv, ...
                         M1p*Lambda*d2qdt2(ret));
dLTdqMGdq  = dLambdaTdq (Lambda,Leq,dLdq,gdofs,ret,slv, ...
                         M1p*Gamma*dqdt(ret));
dLam       = dLambdadq (Lambda,Leq,dLdq,gdofs,ret,slv,d2qdt2(ret));

% ==================== Warning =====================
% Set csflag = 0 if using complex step when calling
% buildMRQLin. Otherwise, csflag = 1 provides
% better performance within dGammadq.
csflag = 0;
[dGam,dGamd] = dGammadq (csflag,q,P,dqdt,ret,slv,gdofs,idofs,idofm,Tb_h, ...
                         Leq,dLdq,Lambda,Gamma,dqdt(ret));

% Partition matrices in the dimension of the partial derivatives.
[dLamp,rr,cr]     = partitionMatrix (dLam,[],slv);
[dGamp,rr,cr]     = partitionMatrix (dGam,[],slv);
[dGamdp,rr,cr]    = partitionMatrix (dGamd,[],slv);

dM   = (LamT*(M1p*dLamp + dMLp) + dLTdqMLd2q)*Lambda;
dMG  = (LamT*(M1p*dGamp + dMGp) + dLTdqMGdq)*Lambda ...
     + LamT*M1p*dGamdp*Gamma;
dMGd = LamT*M1p*(Gamma + dGamdp*Lambda);

% ----------------------------------------------------------------
% RHS vector.
Rvecp = -(G1p*dqdtp - H1p + C1p*dqdtp + K1p);
R = LamT*Rvecp;

dRmat  = -(dG1p - dH1p + dC1p + dK1p);
dRdmat = -(G1p + dGd1p - dHd1p + C1p);

dLTdqR  = dLambdaTdq (Lambda,Leq,dLdq,gdofs,ret,slv,Rvecp);
dR = dLTdqR*Lambda + LamT*(dRmat*Lambda + dRdmat*Gamma);
dRd = LamT*dRdmat*Lambda;











