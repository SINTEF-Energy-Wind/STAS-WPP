function [M,MG,R,Q,Lambda,Gamma,Leq,dLdq] = ...
                             buildMRQNL (s,q,dqdt,P,F)
%
% Assemble the nonlinear structural matrices and vectors.  
%
% Version:        Changes:
% --------        -------------
% 03.02.2017      Original code.
% 18.05.2018      Updated with RHS as output, rather than separate
%                 terms, in order to facilitate matching linear and
%                 nonlinear functions.
% 24.02.2019      Deleted call to hydrodynamic added mass function;
%                 now this is handled as extra mass in the element
%                 mass matrices.
%
% Version:        Verification:
% --------        -------------
% 03.02.2017      
% 18.05.2018      
% 24.02.2019      
%
% Inputs:
% -------
% s               : Data structure defining the turbine structures.
% q,dq/dt         : Present nodal DOFs and their rate of change.
% P               : Undeformed nodal coordinates.
% F               : Nodal forces, including 6 joint DOFs.
%
% Outputs:
% --------
% M               : Matrix that multiplies reduced d2qh/dt2 in the
%                   equations of motion.
% MG              : Lambda^T M Gamma dq/dt term.
% R               : = G-H+C+K, constrained via Lambda^T.
% Q               : Vector of constrained generalized forces.
% Lambda, Gamma   : Constraint matrices.
% Leq             : The L matrix, constraint equation Jacobian.
% dLdq            : dL/dq, second partial derivatives of the constraint
%                   equations.

Ndj = size(q,1);
[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);

[Mz,Gzm,Hz,Kz,Czm] = assembleElementsNL (s,idofs,q,dqdt,P);

% Collapse the Gzm matrix to a force vector.  [Think about rewriting
% assembleElementsNL and buildElementsNL with G as a force vector,
% to reduce storage and speed up the calculation.]
Gz = full(Gzm*dqdt);

% Introduce a damping force vector.
Cz = full(Czm*dqdt);

% ----------------------------------------------------------------
% Add soil stiffness and damping to the foundation.  The function
% returns the nodal forces.  Applying the nodal forces requires the
% operation Qh*F.
%
% [For now, this is based on linear soil springs.  It is implemented
% in terms of a spring force, which can easily be made nonlinear at
% a future time.]
Ndf = 6*s.foundation.Nnod;
dofs = idofs(1) + [1:Ndf].';
Fsoil = soilNL (s.foundation.k(1,:).',s.foundation.k(2,:).',                     ...
                s.foundation.k(3,:).',s.foundation.k(4,:).',s.foundation.k(5,:).', ...
                s.foundation.k(6,:).',s.foundation.k(7,:).',s.foundation.k(8,:).', ...
                P(dofs),q(dofs),dqdt(dofs));

% Add the effect of soil and added mass forces.
Fnod = F;
Fnod(dofs) = Fnod(dofs) + Fsoil;

% ----------------------------------------------------------------
% Apply forces.
Qh = Qhat (q,P,idofs,inods);
QF = Qh*Fnod;

% Compute constraint matrices Lambda, containing the partitioned
% constraint equations, and Gamma, a term involved in computing the
% constraint forces.
[Lambda,Leq,Con,ret,slv] = constraints (q,P,Tb_h,idofs,idofm);
[Gamma,dLdq] = constraintGamma (q,dqdt,P,Leq,Lambda, ...
                                Tb_h,ret,slv,idofs,idofm);
LamT = Lambda.';

% Partition and constrain.
[Mp,rr,cr]  = partitionMatrix (Mz,slv,slv);
[Gp,rr,cr]  = partitionMatrix (Gz,slv,[]);
[Hp,rr,cr]  = partitionMatrix (Hz,slv,[]);
[Cp,rr,cr]  = partitionMatrix (Cz,slv,[]);
[Kp,rr,cr]  = partitionMatrix (Kz,slv,[]);
[Qp,rr,cr]  = partitionMatrix (QF,slv,[]);

M  =  LamT*Mp*Lambda;
MG =  LamT*Mp*Gamma*dqdt(ret);
R  = -LamT*(Gp - Hp + Cp + Kp);
Q  =  LamT*Qp;

%{

'NL'
%full([[1:Ndj].' imag(Qh(90,:)).'.*real(Fnod)/sqrt(eps) ...
%                real(Qh(90,:)).'.*imag(Fnod)/sqrt(eps)])
full([[1:Ndj].' imag(QF)/sqrt(eps)])

fid = fopen('rep.txt','a');
fprintf(fid,'\n');

fprintf(fid,'Q NL\n');
dLQ = imag(LamT)*real(Qp);
LdQ = real(LamT)*imag(Qp);
for ii = 1:6
   fprintf(fid,'%+5.6e %+5.6e %+5.6e\n', ...
           dLQ(72+ii)/sqrt(eps),LdQ(72+ii)/sqrt(eps),imag(Q(72+ii))/sqrt(eps));
end
fprintf(fid,'\n');

fclose(fid);
%}



