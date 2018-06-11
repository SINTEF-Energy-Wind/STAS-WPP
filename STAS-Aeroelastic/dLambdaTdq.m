function dL = dLambdaTdq (Lambda,L,dLdq,dofs,ret,slv,v)
%
% Compute the derivative of the transposed, partitioned/inverted
% constraint equation matrix, times some vector, (dLambda^T/dq)*v.
%
% Version:        Changes:
% --------        -------------
% 18.05.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 18.05.2018      Verified that this matches Lambda using complex
%                 step derivatives and constraints.m.
%
% Inputs:
% -------
% Lambda          : Partitioned/inverted constraint matrix.
% L               : Constraint equations.
% dLdq            : Derivatives of constraint equations.
% dofs            : Master and slave dofs, wrt which dLdq is nonzero.
% ret,slv         : Lists of retained and slave DOFs.
% v               : Vector to postmultiply by.
%
% Outputs:
% --------
% dL              : (dLambda^T/dq)*v

Ndj  = size(Lambda,1);
Nret = size(ret,1);
Nslv = size(slv,1);
Ndof = size(dofs,1);

rtsl = [ret;slv];

dL = spalloc (Nret,Ndj,5*Ndj);
Ls    = L(:,slv);
ILsT   = [speye(Nret) sparse(Nret,Nslv);sparse(Nslv,Nret) Ls].';
ILsTiv = ILsT\v;

for idof = 1:Ndof

   icN = Ndj*(idof-1);


   dLhdq = [sparse(Nret,Nret);dLdq(:,icN+ret)];
   dLsdq = [sparse(Nret,Nret) sparse(Nret,Nslv); ...
            sparse(Nslv,Nret) dLdq(:,icN+slv)];

   dL(:,rtsl==dofs(idof)) = -((dLhdq + dLsdq*Lambda).')*ILsTiv;

end

 