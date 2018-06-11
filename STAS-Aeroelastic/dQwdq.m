function dQ = dQwdq (qn,qB,Pn,PB)
%
% Computes the derivatives of the linear velocity associated with a
% single node.
%
% w = Qw dq/dt
%
% Version:        Changes:
% --------        -------------
% 20.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 20.12.2017      Verified by complex step using Qwnod.m
%
% Inputs:
% -------
% qn              : DOFs associated with the node.
% qB              : DOFs associated with the body reference node.
% Pn              : Undeformed pos,rot of the node wrt body reference.
% PB              : Undeformed pos,rot of the body reference.
%
% Outputs:
% --------
% dQ              : 3-by-12^2 matrix of global velocities, in global
%                   coordinates.  Qw is a function of 12 state
%                   variables: the six body reference DOFs, and the six
%                   nodal DOFs.  So Qw is 3-by-12, and there is a 
%                   partial derivative with respect to each of the 12
%                   DOFs.

dQ = zeros(3,12*12);

r = qn(1:3) + Pn(1:3);

[Fr,TBB0,dTBB0,TB0g] = FRefMatrix (qB(4:6),PB(4:6));
d2TBB0 = d2Tdth2 (qB(4:6),TBB0,dTBB0);

% dQ/dOj = 0.  Translation of the body does not influence any of the
% relationships between v and dO/dt, dTH/dt, dd/dt, or dth/dt.

jref = 3;        % dQ/dTHj
for jdof = 1:3

   jc3 = 3*(jdof-1);
   jc9 = 9*(jdof-1);
   jc12 = 12*(jref+jdof-1);

   dQ(:,jc12+[7:9])   = TB0g*dTBB0(:,jc3+[1:3]);
   dQ(:,jc12+4)       = TB0g*d2TBB0(:,jc9+[1:3])*r;
   dQ(:,jc12+5)       = TB0g*d2TBB0(:,jc9+[4:6])*r;
   dQ(:,jc12+6)       = TB0g*d2TBB0(:,jc9+[7:9])*r;

end

jref = 6;       % dQ/ddn
T1 = TB0g*dTBB0(:,1:3);
T2 = TB0g*dTBB0(:,4:6);
T3 = TB0g*dTBB0(:,7:9);
for jdof = 1:3

   jc12 = 12*(jref+jdof-1);

   dQ(:,jc12+4) = T1(:,jdof);
   dQ(:,jc12+5) = T2(:,jdof);
   dQ(:,jc12+6) = T3(:,jdof);

end

