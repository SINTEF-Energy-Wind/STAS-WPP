function [M,dML,dMG,G,dG,dGd,H,dH,dHd,K,dK] =                 ...
              buildElementLin (mes,kes,dqdt,d2qL,d2qG,        ...
                               Qu1,Qu2,dQu1,dQu2,d2Qu1,d2Qu2, ...
                               mu,dmu,d2mu,TsB,dTsB,d2TsB)
%
% Build the mass, gyroscopic, and stiffness matrices for one element.
% These are each 18-by-18 matrices, as the motion of each element depends
% on its two connected nodes, plus the body reference node.
%
% In general the matrices depend on the pose of the body and its first
% and second time derivatives.  The pose q is implicit in the Qu and TsB
% values input.  The rate and acceleration appear explicitly in the
% equations and must be input.  
%
% The code is developed such that linearized matrices can be obtained
% even around transitory dynamic states (what I'm calling "tangent
% dynamics"), where the "mean" velocity and acceleration are not
% necessarily zero.
%
% Note, the term (dQh/dq) F0 that acts as part of the stiffness matrix
% is not included here at the element level.  It must be incorporated
% later at the structural level.  Damping is also not included here,
% a modal implementation being preferred.
%
% Version:        Changes:
% --------        -------------
% 13.11.2017      Original code.
% 18.05.2018      Updated such that the function is a direct linearization
%                 of the outputs of buildElementNL, rather than combining
%                 the terms into M, C, and K at this stage.
% 15.04.2019      Restructured the calculations to accelerate the for loops.
%
% Version:        Verification:
% --------        -------------
% 13.11.2017    - The derivatives of M, G, K with respect to q, qdot
%                 have all been verified by complex step, using a
%                 nonlinear version of the equations, beamElementNL.m.
%               - This function has been used to successfully solve the
%                 deflection of a 45-degree circular bend subject to an
%                 out-of-plane load.  For instance, Bathe and Bolourchi
%                 (1979) p 982.  This validates the approach to TsB and
%                 mu.
% 18.05.2018      Derivatives verified again with complex step using
%                 buildElementNL.  [CHECK VELOCITY TERMS]
% 15.04.2019      
%
% Inputs:
% -------
% mes, kes        : Element mass and stiffness matrices.
% dqdt            : Rate of change of DOFs.
% d2qL            : Acceleration term with constraint matrix, Lambda*d2qdt2.
% d2qG            : Gamma*dqdt.
% Qu1,2           : [v;w] = Qu*q.  Qu is 6-by-12.  The latter dimension is
%                   6 body reference node and 6 nodal DOFs.
% dQu1,2          : A 6-by-12*12 matrix containing dQu/dq for each node.
% d2Qu1,2         : A 6-by-12*12*12 matrix of second derivatives.
% mu,dmu,d2mu     : Elastic deformations multiplying the element stiffness
%                   matrix -- and their derivatives.
% TsB,dTsB,d2TsB  : Section-to-body transform and its derivatives.
%
% Outputs:
% --------
% M               : 18-by-18 mass matrix: qB, qn1, qn2 DOFs.
% dM              : 18-by-18 matrix: dM/dq*d2q0/dt2.
% G               : 18-by-18 gyroscopic matrix.
% dG              : 18-by-18 matrix: dG/dq*dq0/dt.
% dGd             : 18-by-18 matrix: dG/dqdot*dq0/dt.
% H               : 18-by-1 centrifugal vector.
% dH              : 18-by-18 matrix: dH/dq.
% dHd             : 18-by-18 matrix: dH/dqdot.
% K               : 18-by-1 elastic stiffness vector.
% dK              : 18-by-18 matrix: dK/dq.

% 12-by-12, 12-by-12*18, 12-by-12*18*18.
[TmT,dTmT,d2TmT] = d2TmeT (mes,TsB,dTsB,d2TsB);

% Form the two nodal Qu matrices into a format that multiplies a 12-by-12
% element matrix.   Que is 12-by-18.
%               |O  |
%  |v1|         |Phi|
%  |w1| = [Que] |d1 |
%  |v2|         |th1|
%  |w2|         |d2 |
%               |th2|
Que   = Qel (Qu1,Qu2);        % 12-by-18.
dQue  = dQel (dQu1,dQu2);     % 12-by-18*18.
d2Que = d2Qel (d2Qu1,d2Qu2);  % 12-by-18*18*18.

Qp = Que.';
dQp = dQue.';
d2Qp = d2Que.';

% Terms to be computed upfront.
ue = Que*dqdt;
uep = ue.';
dqdtp = dqdt.';

TmTQ = TmT*Que;
TmTu = TmT*ue;

dQTmTQ = dQp*TmTQ;

% --------------------------------------------------------------------
% Mass matrix.
M = sparse((Que.')*TmTQ);

% Other terms with the mass matrix.  Accelerated and verified.
dML = zeros(18,18);
dMG = zeros(18,18);
mat = zeros(18,18*18);
for kk = 1:18

   kc18 = 18*(kk-1);
   kc12 = 12*(kk-1);

   mat(:,kc18+[1:18]) = dQTmTQ(kc18+[1:18],:) + dQTmTQ(kc18+[1:18],:).' ...
                      + (Que.')*dTmT(:,kc12+[1:12])*Que;
   dML(:,kk) = mat(:,kc18+[1:18])*d2qL;
   dMG(:,kk) = mat(:,kc18+[1:18])*d2qG;

end

% --------------------------------------------------------------------
% Gyroscopic matrix G.  Accelerated and verified.
G = zeros(18,18);
for jj = 1:18

   jc18 = 18*(jj-1);
   jc12 = 12*(jj-1);

   G(:,jj) = mat(:,jc18+[1:18])*dqdt;

end

% Accelerated and verified.
dG = zeros(18,18);
dQpdTmT = dQp*dTmT;
dQpTmTdQ = dQp*TmT*dQue;
d2QpTmTQ = d2Qp*TmT*Que;
d2QpTmTu = d2Qp*TmTu;
for jj = 1:18

   jc324 = 324*(jj-1);
   jc216 = 216*(jj-1);
   jc18  =  18*(jj-1);
   jc12  =  12*(jj-1);

   for kk = 1:18

      kc18 = 18*(kk-1);
      kc12 = 12*(kk-1);

      dG(:,kk) = dG(:,kk)                                              ...
               + (d2QpTmTu(jc324+kc18+[1:18])                          ...
               +  dQpdTmT(jc18+[1:18],kc12+[1:12])*ue                  ...
               +  dQpTmTdQ(jc18+[1:18],kc18+[1:18])*dqdt               ...
               +  (d2QpTmTQ(jc324+kc18+[1:18],:).')*dqdt               ...
               +  dQpTmTdQ(kc18+[1:18],jc18+[1:18])*dqdt               ...
               +  Qp*(dQpdTmT(jc18+[1:18],kc12+[1:12]).')*dqdt         ...
               +  dQpdTmT(kc18+[1:18],jc12+[1:12])*ue                  ...
               +  Qp*(dQpdTmT(kc18+[1:18],jc12+[1:12]).')*dqdt         ...
               +  Qp*d2TmT(:,jc216+kc12+[1:12])*ue)*dqdt(jj);

   end

end

% Accelerated and verified.
dGd = zeros(18,18);
for jj = 1:18

   jc18 = 18*(jj-1);
   jc12 = 12*(jj-1);

   dGd = dGd + mat(:,jc18+[1:18])*dqdt(jj);

end

% --------------------------------------------------------------------
% Centrifugal vector H.
H = zeros(18,1);
for jj = 1:18

   jc18 = 18*(jj-1);
   jc12 = 12*(jj-1);

   H(jj) = H(jj)                               ...
         + dqdtp*dQp(jc18+[1:18],:)*TmTu       ...
         + 0.5*uep*dTmT(:,jc12+[1:12])*ue;

end

% Accelerated and verified.
dH = zeros(18,18);
for ii = 1:18

   ic324 = 324*(ii-1);
   ic216 = 216*(ii-1);
   ic18  =  18*(ii-1);
   ic12  =  12*(ii-1);

   for kk = 1:18

      kc18 = 18*(kk-1);
      kc12 = 12*(kk-1);

      % Verified.
      dH(ii,kk) = dH(ii,kk)                                     ...
                + dqdtp*d2QpTmTu(ic324+kc18+[1:18])             ...
                + dqdtp*dQpTmTdQ(ic18+[1:18],kc18+[1:18])*dqdt  ...
                + dqdtp*dQpdTmT(ic18+[1:18],kc12+[1:12])*ue     ...
                + dqdtp*dQpdTmT(kc18+[1:18],ic12+[1:12])*ue     ...
                + 0.5*uep*d2TmT(:,ic216+kc12+[1:12])*ue;

   end

end

% Accelerated and verified.
dHd = zeros(18,18);
for ii = 1:18

   ic18 = 18*(ii-1);

   % Verified.
   dHd(ii,:) = dHd(ii,:) + dqdtp*mat(:,ic18+[1:18]);

end

% --------------------------------------------------------------------
% Stiffness vector K.
K = zeros(18,1);
for jj = 1:12
   K(jj+6) = K(jj+6) + (dmu(:,jj).')*kes*mu;
end

% Accelerated and verified.
dK = zeros(18,18);
dK(7:18,7:18) = (dmu.')*kes*dmu;
kesmu = kes*mu;
for ii = 1:12

   ic12 = 12*(ii-1);

   dK(7:18,ii+6) = dK(7:18,ii+6) + (d2mu(:,ic12+[1:12]).')*kesmu;

end

% (Note that dG is antisymmetric, which destroys the symmetry
% properties of the final assembled K matrix.)


