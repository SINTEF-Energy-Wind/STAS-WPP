function Q = Qwnod (qn,qB,Pn,PB)
%
% Computes the linear velocity, associated with a single node, in
% global coordinates.
%
% w_n/g^g = Qw*dq/dt.  This function builds Qw.
%
% Version:        Changes:
% --------        -------------
% 20.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 20.12.2017      Checked some cases.
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

Q = zeros(3,12);

r = qn(1:3) + Pn(1:3);

[Fr,TBB0,dTBB0,TB0g] = FRefMatrix (qB(4:6),PB(4:6));
TBg = TB0g*TBB0;

Q(:,7:9)   = TBg;
Q(:,1:3)   = eye(3);
Q(:,4)     = TB0g*dTBB0(:,1:3)*r;
Q(:,5)     = TB0g*dTBB0(:,4:6)*r;
Q(:,6)     = TB0g*dTBB0(:,7:9)*r;

