function [LHS,A,Bu,By,C,Du,Dy,                         ...
          Lambda,Gamma,shape,freq,mdamp,ret,slv] =     ...
               structureLin (bmopt,shpflag,            ...
                             s,q0,dq0dt,d2qindt2,P,F0, ...
                             shapein,mdampin,grav)
%
% Assemble the structural equations in state space form.  Note that
% q = (Lamq shape) eta, where Lamq has undone the partitioning of 
% Lambda.
%
%   States:              y vector:             u vector:
%   eta        N         q          Ndj        F          Ndj
%   deta/dt    N         dqdt       Ndj
%                        d2qdt2     Ndj
%
% Version:        Changes:
% --------        -------------
% 22.05.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 22.05.2018      Match with structureNL verified by complex step.
%                 Converged Newton solutions for nonlinear equations.
%                 NREL 5 MW turbine modes match published values.
%
% Inputs:
% -------
% bmopt           : 1 to apply body mode reduction, 0 to use the full
%                   matrices.
% shpflag         : 1 to compute shape, damping etc from the present
%                   matrices, 0 to use existing basis functions, via
%                   shapein and mdampin.
% s               : Data structure describing the structure.
% q0...d2q0dt2    : DOF and first two time derivatives.
% P               : Undeformed nodal positions.
% F0              : Ndj nodal forces.
% shapein,mdampin : Mode shapes and modal damping vector to use if
%                   shpflag = 0.
% grav            : Gravity vector in global coordinates.
%
% Outputs:
% --------
% A...Dy          : State matrices.
% Lambda, Gamma   : Constraint matrices.
% shape,freq      : Body modes and frequencies.
% ret,slv         : Retained and slave DOFs.

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

Ndj = size(q0,1);

gvec = zeros(Ndj,1);
gvec(1:3) = grav;
d2q0dt2 = d2qindt2 - gvec;

if (bmopt == 1)
   % Use the body modes as a basis for DOF reduction.  In these
   % outputs, dL2 = (dLam/dq)*(d2qh0/dt2) and dG2 = (dGam/dq)*(dqh0/dt).
   % The rows are partitioned according to [ret;slv] as is the
   % standard for Lambda and Gamma.  The columns, corresponding to the
   % derivative variable q, are not partitioned, they have the
   % original ordering of the Ndj dofs.
   [M,D,K,slv,ret,Lambda,Gamma,                               ...
    Leq,dLeqdq,dLam2,dGam1,dGamd1,shape,freq,mdamp] =         ...
                        buildMCK (shpflag,s,q0,dq0dt,d2q0dt2, ...
                                  P,F0,shapein,mdampin);

elseif (bmopt == 0)
   % Use the full matrices.
   [M1,dM1,MG1,dMG1,dMGd1,R1,dR1,dRd1,Q1,dQ1,dQd1,slv,ret, ...
    Lambda,Gamma,Leq,dLeqdq,dLam2,dGam1,dGamd1] =          ...
                        buildMRQLin (s,q0,dq0dt,d2q0dt2,P,F0);
   shape = speye(size(M1,1));
   freq  = zeros(size(M1,1),1);
   mdamp = zeros(size(M1,1),1);

   M = M1;
   D = dMGd1 - dRd1 - dQd1;
   K = dM1 + dMG1 - dR1 - dQ1; 

end

N = size(M,1);

% Put these into state space.  This first definition can be thought
% of as building a state space of the modal-transformed equations.
%LHS = [speye(N,N) sparse(N,N);sparse(N,N) M];
%A   = [sparse(N,N) speye(N,N);-K -D];
%
% This normalization is the result of doing a modal transformation
% AFTER building the state space.
sTs = (shape.')*shape;
LHS = [sTs sparse(N,N);sparse(N,N) M];
A   = [sparse(N,N) sTs;-K -D];

LA = LHS\A;  % For nodal accelerations, but watch out for conditioning.

rsind         = [ret;slv];
gdofs         = getgdofs (idofs,idofm,Ndj);
Lamq          = sparse (size(Lambda,1),size(Lambda,2));
Lamq(rsind,:) = Lambda;         % Undo the partitioning, so q = Lamq*qhat.
dLq           = sparse (size(dLam2,1),size(dLam2,2));
dLq(rsind,:)  = dLam2;
LSh           = Lamq*shape;
LShT          = LSh.';
dLam1         = sparse (Ndj,Ndj);
dLam1(rsind,:)= dLambdadq (Lambda,Leq,dLeqdq,gdofs,ret,slv,dq0dt(ret));
Gamq          = sparse (size(Gamma,1),size(Gamma,2));
Gamq(rsind,:) = Gamma;          % d2q/dt2 = Lamq*d2qh/dt2 + Gamq*dqh/dt.
dGq           = sparse (size(dGam1,1),size(dGam1,2));
dGq(rsind,:)  = dGam1;
dGdq          = sparse (size(dGamd1,1),size(dGamd1,2));
dGdq(rsind,:) = dGamd1;
GSh           = Gamq*shape;

Qh = Qhat (q0,P,idofs,inods);
By = sparse(2*N,3*Ndj);
Bu = [sparse(N,Ndj);LShT*Qh];
LBu = LHS\Bu;

% We have the following linearized equations:
% dq = Lamq*shape*deta.  dq/dt = Lamq*shape*deta/dt + Gamq*shape*deta.
% d2q/dt2 = Lamq*shape*d2eta/dt2 
%         + Gamq*shape*deta/dt  + dGam/dqdot*deta0/dt*dq/dt
%         + (dLam/dq*d2eta0/dt2 + dGam/dq*deta0/dt)*dq.
% d2eta/dt2 = (M^-1)(-K deta - D deta/dt + LShT*Qh*F)
%
% So we have LSh and GSh with rows in Ndj order.  The rows and
% columns of dLq and dGq are also in Ndj order.
C  = [LSh            sparse(Ndj,N); ...
      GSh            LSh;           ...
      LSh*LA(N+[1:N],1:N) LSh*LA(N+[1:N],N+[1:N])+GSh];
Du = [sparse(2*Ndj,Ndj);LSh*LBu(N+[1:N],:)];
Dy = [sparse(Ndj,3*Ndj); ...
      sparse(Ndj,3*Ndj); ...
      dLq+dGq dGdq sparse(Ndj,Ndj)];

