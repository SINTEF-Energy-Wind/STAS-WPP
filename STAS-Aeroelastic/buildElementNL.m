function [M,G,H,K] = buildElementNL (mes,kes,dqdt,      ...
                                     Qu1,Qu2,dQu1,dQu2, ...
                                     mu,dmu,TsB,dTsB)
%
% Build the nonlinear mass matrix M, gyroscopic matrix G, centrifugal
% force vector H, and stiffness force vector K, for a single element.
% The matrices are each 18-by-18, and the vectors are of length 18, as
% the motion of each element depends on its two connected nodes, plus
% the body reference node.
%
% In general the matrices depend on the pose of the body and its first
% and second time derivatives.  The pose q is implicit in the Qu and TsB
% values input.  The rate and acceleration appear explicitly in the
% equations and must be input.
%
% Version:        Changes:
% --------        -------------
% 23.11.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 23.11.2017      Complex step derivatives of these equations match the
%                 linearized terms in buildElementLin.m.
%
% Inputs:
% -------
% mes, kes        : Element mass and stiffness matrices.
% dqdt            : Rate of change of DOFs.
% Qu1,2           : [v;w] = Qu*q.  Qu is 6-by-12.  The latter dimension is
%                   6 body reference node and 6 nodal DOFs.
% dQu1,2          : A 6-by-12*12 matrix containing dQu/dq for each node.
% mu,dmu          : Elastic deformations multiplying the element stiffness
%                   matrix -- and their derivatives.
% TsB,dTsB        : Section-to-body transform and its derivatives.
%
% Outputs:
% --------
% M               : 18-by-18 mass matrix: qB, qn1, qn2 DOFs.
% G               : 18-by-18 gyroscopic matrix.
% H               : 18-by-1 centrifugal vector.
% K               : 18-by-1 elastic stiffness vector.

%'buildElementNL'

% 12-by-12, 12-by-12*18.
[TmT,dTmT] = dTmeT (mes,TsB,dTsB);

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

% Velocity terms to be computed upfront.
ue = Que*dqdt;
uep = ue.';

TmTQ = TmT*Que;
TmTu = TmT*ue;

% --------------------------------------------------------------------
% Mass matrix. 
M = zeros(18,18);
M = (Que.')*TmT*Que;

% --------------------------------------------------------------------
% Gyroscopic matrix G.
G = zeros(18,18);
for jj = 1:18

   jc18 = 18*(jj-1);
   jc12 = 12*(jj-1);

   G(:,jj) = G(:,jj)                            ...
           + (dQue(:,jc18+[1:18]).')*TmTu       ...
           + (TmTQ.')*dQue(:,jc18+[1:18])*dqdt  ...
           + (Que.')*dTmT(:,jc12+[1:12])*ue;

end

% --------------------------------------------------------------------
% Centrifugal vector H.
H = zeros(18,1);
for jj = 1:18

   jc18 = 18*(jj-1);
   jc12 = 12*(jj-1);

   H(jj) = H(jj)                               ...
         + ((dQue(:,jc18+[1:18])*dqdt).')*TmTu ...
         + 0.5*uep*dTmT(:,jc12+[1:12])*ue;

end

% --------------------------------------------------------------------
% Stiffness vector K.
K = zeros(18,1);
for jj = 1:12
   K(jj+6) = K(jj+6) + (dmu(:,jj).')*kes*mu;
end

M = sparse(M);
G = sparse(G);

