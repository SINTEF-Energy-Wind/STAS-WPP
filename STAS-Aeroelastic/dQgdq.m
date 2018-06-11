function dQ = dQgdq (qn,qB,Pn,PB)
%
% Computes the derivatives of the linear and angular velocities
% associated with a single node.
%
% [v;w] = Qu dq/dt
%
% Version:        Changes:
% --------        -------------
% 20.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 20.12.2017      Verified by complex step on Qgnod.m.
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
% dQ              : Global velocity, in body coordinates.  6 rows
%                   representing [v;w].  Qu is a function of 12 state
%                   variables: the six body reference DOFs, and the six
%                   nodal DOFs.  So Qu is 6-by-12, and there is a 
%                   partial derivative with respect to each of the 12
%                   DOFs.

dQ = zeros(6,12*12);

r = qn(1:3) + Pn(1:3);

[Fr,TBB0,dTBB0,TB0g] = FRefMatrix (qB(4:6),PB(4:6));
d2TBB0 = d2Tdth2 (qB(4:6),TBB0,dTBB0);
Fmat = Dmatrix (Fr);

TBg = TB0g*TBB0;
Tn0B = TFromTheta (Pn(4:6));

[F,Tnn0,dTnn0] = Fmatrix (qn(4:6),TBg,Tn0B);
d2Tnn0 = d2Tdth2 (qn(4:6),Tnn0,dTnn0);
Gmat = Dmatrix (F);

dFr = dFrefdth (qB(4:6),TBB0,dTBB0,d2TBB0,TB0g);
dFmat = zeros(3,3*3);
for jj = 1:3
   jc3 = 3*(jj-1);
   jc9 = 9*(jj-1);
   dFmat(:,jc3+[1:3]) = Dmatrix (dFr(:,jc9+[1:9]));
end

dF = dFdth (Tnn0,dTnn0,d2Tnn0,TBB0,dTBB0,TB0g,Tn0B);
dGmat = zeros(3,6*3);
for jj = 1:6
   jc3 = 3*(jj-1);
   jc9 = 9*(jj-1);
   dGmat(:,jc3+[1:3]) = Dmatrix (dF(:,jc9+[1:9]));   
end

% dQ/dOj = 0.  Translation of the body does not influence any of the
% relationships between [v;w] and dO/dt, dTH/dt, dd/dt, or dth/dt.

jref = 3;        % dQ/dTHj
for jdof = 1:3

   jc3 = 3*(jdof-1);
   jc9 = 9*(jdof-1);
   jc12 = 12*(jref+jdof-1);

   dQ(1:3,jc12+[7:9])   = TB0g*dTBB0(:,jc3+[1:3]);
   dQ(1:3,jc12+4)       = TB0g*(d2TBB0(:,jc9+[1:3])*r);
   dQ(1:3,jc12+5)       = TB0g*(d2TBB0(:,jc9+[4:6])*r);
   dQ(1:3,jc12+6)       = TB0g*(d2TBB0(:,jc9+[7:9])*r);
   dQ(4:6,jc12+[4:6])   = dFmat(:,jc3+[1:3]);
   dQ(4:6,jc12+[10:12]) = dGmat(:,jc3+[1:3]);

end

jref = 6;       % dQ/ddn
T1 = TB0g*dTBB0(:,1:3);
T2 = TB0g*dTBB0(:,4:6);
T3 = TB0g*dTBB0(:,7:9);
for jdof = 1:3

   jc12 = 12*(jref+jdof-1);

   dQ(1:3,jc12+4) = T1(:,jdof);
   dQ(1:3,jc12+5) = T2(:,jdof);
   dQ(1:3,jc12+6) = T3(:,jdof);

end

jref = 9;       % dQ/dthn
for jdof = 1:3

   jc3 = 3*(jdof-1);
   jc12 = 12*(jref+jdof-1);

   dQ(4:6,jc12+[10:12]) = dGmat(:,9+jc3+[1:3]);

end
