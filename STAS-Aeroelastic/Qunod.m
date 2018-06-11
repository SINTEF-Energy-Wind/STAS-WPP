function Q = Qunod (qn,qB,Pn,PB)
%
% Computes the linear and angular velocities associated with a single
% node.
%
% [v^B;w^B] = Qu*dq/dt.  This function builds Qu.
%
% Version:        Changes:
% --------        -------------
% 26.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 26.10.2017      Checked in-depth on a handful of cases.
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
% Q               : Matrix giving global velocity, in body coordinates.

Q = zeros(6,12);

r = qn(1:3) + Pn(1:3);

[Fr,TBB0,dTBB0,TB0g] = FRefMatrix (qB(4:6),PB(4:6));
Fmat = Dmatrix (Fr);
TBg = TB0g*TBB0;
TgB = TBg.';
Tn0B = TFromTheta (Pn(4:6));
[F,Tnn0,dTnn0] = Fmatrix (qn(4:6),TBg,Tn0B);
Gmat = Dmatrix (F);

Q(1:3,7:9)   = eye(3);
Q(1:3,1:3)   = TgB;
Q(1:3,4)     = TgB*TB0g*dTBB0(:,1:3)*r;
Q(1:3,5)     = TgB*TB0g*dTBB0(:,4:6)*r;
Q(1:3,6)     = TgB*TB0g*dTBB0(:,7:9)*r;
Q(4:6,4:6)   = TgB*Fmat;
Q(4:6,10:12) = TgB*Gmat;
