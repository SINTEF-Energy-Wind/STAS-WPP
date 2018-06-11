function C = nodalDOFs (shape,Lambda,ret,slv)
%CHECK THIS  Q = C ETA IS INCORRECT, I THINK.
%
% It is often convenient to have the unconstrained, un-reduced nodal DOFs.
% These are back-calculated from the states by undoing the reduction,
% constraint, and partitioning procedures.
%
% q = C eta,  dq/dt = C deta/dt.
%
%   States:           y vector:         u vector:
%   -------           ---------         ---------
%   eta               q
%   deta/dt           dq/dt
%
% Version:        Changes:
% --------        -------------
% 21.01.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 21.01.2017      
%
% Inputs:
% -------
% shape           : Mode shapes, Nrdof-by-Nx.
% Lambda          : Constraint matrix, Ndof-by-Nrdof, Ndof partitioned into
%                   [ret;slv].
% ret,slv         : Retained and slave DOFs, out of Ndof total.
%
% Outputs:
% --------
% C               : State matrices.

Nx    = size (shape,2);
Nrdof = size (ret,1);
Ndof  = size (Lambda,1);

C = spalloc (2*Ndof,2*Nx,2*(Nx*Ndof));

C(ret,1:Nx) = Lambda(1:Nrdof,:)*shape;
C(slv,1:Nx) = Lambda(Nrdof+1:Ndof,:)*shape;
C(Ndof+[1:Ndof],Nx+[1:Nx]) = C(1:Ndof,1:Nx);

