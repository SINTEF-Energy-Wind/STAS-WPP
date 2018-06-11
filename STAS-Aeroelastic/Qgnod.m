function Q = Qgnod (qn,qB,Pn,PB)
%
% Computes the linear and angular velocities associated with a single
% node.
%
% [v^g;w^g] = Qu*dq/dt.  This function builds Qu.
%
% Version:        Changes:
% --------        -------------
% 20.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 20.12.2017      Checked in-depth on a handful of cases.
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
% Q               : Matrix giving global velocity, in global coordinates.

Q = zeros(6,12);

r = qn(1:3) + Pn(1:3);

[Fr,TBB0,dTBB0,TB0g] = FRefMatrix (qB(4:6),PB(4:6));
Fmat = Dmatrix (Fr);
TBg = TB0g*TBB0;
Tn0B = TFromTheta (Pn(4:6));
[F,Tnn0,dTnn0] = Fmatrix (qn(4:6),TBg,Tn0B);
Gmat = Dmatrix (F);

Q(1:3,1:3)   = eye(3);
Q(1:3,7:9)   = TBg;
Q(1:3,4)     = TB0g*dTBB0(:,1:3)*r;
Q(1:3,5)     = TB0g*dTBB0(:,4:6)*r;
Q(1:3,6)     = TB0g*dTBB0(:,7:9)*r;
Q(4:6,4:6)   = Fmat;
Q(4:6,10:12) = Gmat;

%'---Qgnod---'
%[qn qB Pn PB]
%Q