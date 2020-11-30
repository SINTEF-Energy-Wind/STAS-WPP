function [q,qd,qdd,T,dTdq,dTdqd] = buildTqhq (linFlag,s,qh,qhd,qhdd,P,ret,slv)
%
% Builds the transform that implements the constraints on the
% structural state variables.  Also computes the values of the
% full q vector, solving the nonlinear equations for qs in 
% terms of qh.
%
% d |q |  = T_qh^q d |qh | , with T_qh^q = |Lam  0 | 
% dt|qd|           dt|qhd|                 |Gam Lam|   
%
% Version:        Changes:
% --------        -------------
% 13.01.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 13.01.2020      Verified that q, qd, qdd match the true values,
%                 from which input qh, qhd, qhdd were derived.
%                 Verified T by complex step on qh, qhd inputs.
%
% Inputs:
% -------
% linFlag         : = 1 to include linearization, dTdq and dTdqd.
% s               : Data structure with structural inputs.
% qh,qhd,qhdd     : q^h and its time derivatives.
% P               : Vector of nodal positions.
% ret,slv         : Retained and slave DOFs.
%
% Outputs:
% --------
% q, qd, qdd      : q and its time derivatives.
% T               : T_qh^q.
% dTdq            : (dT/dq)*[qhd;qhdd], need to multiply by the vector
%                   to prevent storing a large number of matrices.
% dTdqd           : (dT/dqd)*[qhd;qhdd]

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
Ndj = size(ret,1) + size(slv,1);

% Get the elements of q that are active in the constraint equations.
gdofs = getgdofs (idofs,idofm,Ndj);

Nret = size(ret,1);   % Number of retained DOFs in each of qh, qhd, qhdd.
II = speye(Nret);

% Will need Tb_h below.
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);

% Call a Newton-Raphson solution for the nonlinear constraint
% equations.  We are not at this point dealing with mode shapes,
% so just set the shape matrix to I.
qguess = zeros(size(slv,1),1);  % Initial guess for Newton-Raphson.
[q,qd] = qFromx ([qh;qhd],s,P,II,ret,slv,qguess);

% Get Lambda and Gamma and supporting variables.  Lambda and Gamma
% are in a partitioned ordering.
[Lambda,Leq,Con,ret,slv] = constraints (q,P,Tb_h,idofs,idofm);
[Gamma,dLdq] = constraintGamma (q,qd,P,Leq,Lambda,Tb_h, ... 
                                ret,slv,idofs,idofm);

% Reorder Lambda and Gamma so partitioning isn't necessary.
rsind         = [ret;slv];
Lamq          = sparse (size(Lambda,1),size(Lambda,2));
Lamq(rsind,:) = Lambda;
Gamq          = sparse (size(Gamma,1),size(Gamma,2));
Gamq(rsind,:) = Gamma;

qdd = Lamq*qhdd + Gamq*qhd;
T = [Lamq, sparse(Ndj,Nret); Gamq, Lamq];

if (linFlag == 1)

   % We are going to need (dLam/dq)*d2qh/dt2, dLam/dq*dqh/dt,
   % dGam/dq*dqh/dt, and dGam/dqd*dqh/dt.
   dLamqhd  = dLambdadq (Lambda,Leq,dLdq,gdofs,ret,slv,qhd);
   dLamqhdd = dLambdadq (Lambda,Leq,dLdq,gdofs,ret,slv,qhdd);

   dLamdq = sparse(Ndj,Ndj);
   dLamdq(rsind,:) = dLamqhd;
   dLamddq = sparse(Ndj,Ndj);
   dLamddq(rsind,:) = dLamqhdd;

   % ==================== Warning =====================
   % Set csflag = 0 if using complex step for gradients.
   % Otherwise, csflag = 1 provides better performance
   % within dGammadq.
   if (isreal(q) && isreal(P) && isreal(qd) && isreal(Tb_h) && ...
       isreal(Leq) && isreal(dLdq) && isreal(Lambda) && ...
       isreal(Gamma) && isreal(qhd))
      csflag = 1;  % No complex step.
   else
      csflag = 0;  % Some imaginary input, could be complex step.
   end
   [dGam,dGamd] = dGammadq (csflag,q,P,qd,ret,slv,gdofs,idofs,idofm,Tb_h, ...
                            Leq,dLdq,Lambda,Gamma,qhd);

   dGamq = sparse(Ndj,Ndj);
   dGamq(rsind,:) = dGam;
   dGamdq = sparse(Ndj,Ndj);
   dGamdq(rsind,:) = dGamd;

   dTdq = [dLamdq; dGamq + dLamddq];
   dTdqd = [sparse(Ndj,Ndj); dGamdq];

else

   dTdq  = sparse (2*Ndj,Ndj);
   dTdqd = sparse (2*Ndj,Ndj);

end


