function dL = dLambdadq (Lambda,L,dLdq,dofs,ret,slv,v)
%
% Compute the derivative of the partitioned/inverted constraint
% equation matrix, times some vector, (dLambda/dq)*v.
%
% Version:        Changes:
% --------        -------------
% 22.03.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 22.03.2018      
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
% dL              : (dLambda/dq)*v

Ndj  = size(Lambda,1);
Nret = size(ret,1);
Nslv = size(slv,1);
Ndof = size(dofs,1);

dL = spalloc (Ndj,Ndj,5*Ndj);
Ls    = L(:,slv);
ILs   = [speye(Nret) sparse(Nret,Nslv);sparse(Nslv,Nret) Ls];

for idof = 1:Ndof

   icN = Ndj*(idof-1);
   
   dLhdq = [sparse(Nret,Nret);dLdq(:,icN+ret)];
   dLsdq = [sparse(Nret,Nret) sparse(Nret,Nslv); ...
            sparse(Nslv,Nret) dLdq(:,icN+slv)];

   dL(:,dofs(idof)) = -ILs\((dLhdq + dLsdq*Lambda)*v);

end

 