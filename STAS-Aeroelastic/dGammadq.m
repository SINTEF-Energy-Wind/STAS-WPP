function [dG,dGd] = dGammadq (csflag,q,P,dqdt,ret,slv,dofs,idofs,idofm,Tbh, ...
                              L,dLdq,Lambda,Gamma,v)
%
% Compute the derivative of the partitioned/inverted constraint
% equation matrix, times some vector, (dGamma/dq)*v.
%
% Version:        Changes:
% --------        -------------
% 22.03.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 22.03.2018      Verified against constraintGamma using complex step.
%
% Inputs:
% -------
% csflag          : = 1 to employ complex step to determine d2L/dq2.
%                   = 0 to employ finite difference.
%                   Use csflag = 0 if dGammadq is being called with
%                   EXTERNAL complex step derivatives.  Otherwise,
%                   csflag = 1 is both faster and more accurate.
% q,P,dqdt        : Nodal quantities.
% ret,slv         : Lists of retained and slave DOFs.
% dofs            : Master and slave dofs, wrt which dLdq is nonzero.
% idofs,idofm     : Slave and master dof indices.
% L               : Constraint equations.
% dLdq            : Derivatives of constraint equations.
% Lambda,Gamma    : Partitioned/inverted constraint matrices.
% v               : Vector to postmultiply by.
%
% Outputs:
% --------
% dG,dGd          : (dGamma/dq)*v,(dGamma/dqdot)*v

Ndj  = size(Lambda,1);
Nret = size(ret,1);
Nslv = size(slv,1);
Ndof = size(dofs,1);

dG  = spalloc (Ndj,Ndj,10*Ndj);
dGd = spalloc (Ndj,Ndj,5*Ndj);

dLamv = dLambdadq (Lambda,L,dLdq,dofs,ret,slv,v);
Lamv  = Lambda*v;
Gamv  = Gamma*v;

rsind = [ret;slv];

for idof = 1:Ndof

   icN = Ndj*(idof-1);

   % Use numerical derivatives for the second derivative of the constraint
   % equations L.  This is difficult to formulate analytically, and
   % does not play a dominant role in the turbine dynamics, so a
   % moderate loss in numerical precision should be acceptable here.
   % (Complex step retains near machine precision, but can't be used if
   % we're already using complex step derivatives from a higher-level
   % calling function.)
   if (csflag == 1)
      del = sqrt(eps);
      qc = q;
      qc(dofs(idof)) = qc(dofs(idof)) + i*del;
      dLdqc = derivConstraints (qc,P,Tbh,idofs,idofm);
      d2L = imag(dLdqc)/del;
   elseif (csflag == 0)
      del = 1e-5;
      qh = q;
      qh(dofs(idof)) = qh(dofs(idof)) + del;
      ql = q;
      ql(dofs(idof)) = ql(dofs(idof)) - del;
      dLdqh = derivConstraints (qh,P,Tbh,idofs,idofm);
      dLdql = derivConstraints (ql,P,Tbh,idofs,idofm);
      d2L   = (dLdqh - dLdql)/(2*del);
   end

   H = zeros (Ndj,1);
   for jdof = 1:Ndof
      jcN = Ndj*(jdof-1);
      H(Nret+[1:Nslv]) = H(Nret+[1:Nslv])                          ...
                - d2L(:,jcN+rsind)*(Lamv*dqdt(dofs(jdof)))         ...
                - dLdq(:,jcN+rsind)*(dLamv(:,dofs(idof))*dqdt(dofs(jdof)));
   end

   H = H - [sparse(Nret,Nret) sparse(Nret,Nslv);     ...
            sparse(Nslv,Nret) dLdq(:,icN+slv)]*Gamv;

   Ls    = L(:,slv);
   ILs   = [speye(Nret) sparse(Nret,Nslv);sparse(Nslv,Nret) Ls];

   dG(:,dofs(idof)) = ILs\H;

   dGd(:,dofs(idof)) = [sparse(Nret,1);-Ls\(dLdq(:,icN+rsind)*Lamv)];

end

 