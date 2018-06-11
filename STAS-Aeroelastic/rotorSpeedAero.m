function [W,azi,Dy] = rotorSpeedAero (q,dqdt,P,Try,Ydof,Ddof,nodof)
%
% Rotor speed and effective azimuth, as experienced by the
% aerodynamic model.
%
%   States:           y vector:         u vector:
%                     q         Ndj
%                     dq/dt     Ndj
%
% Version:        Changes:
% --------        -------------
% 08.02.2018      Original code.
% 03.05.2018      Added calculation of effective azimuth.
%
% Version:        Verification:
% --------        -------------
% 08.02.2018      Linearization verified by complex step.
% 03.05.2018      
%
% Inputs:
% -------
% Ddof            : Driveshaft reference DOF.
% nodof           : A reference DOF for a node on the driveshaft
%                   representing rotor speed.  Say, the hub node.
%
% Outputs:
% --------
% 

Ndj = size(q,1);

dofs = [Ydof+[4:6] Ddof+[4:6] nodof+[4:6]].';
qq   = q(dofs);
PP   = P(dofs);
dqq  = dqdt(dofs);

Dy = spalloc(2,2*Ndj,36);

[Tyy0,dTyy0] = dTdth (qq(1:3));
[Tdd0,dTdd0] = dTdth (qq(4:6));
[Tnn0,dTnn0] = dTdth (qq(7:9));
d2Tyy0 = d2Tdth2 (qq(1:3),Tyy0,dTyy0);
d2Tdd0 = d2Tdth2 (qq(4:6),Tdd0,dTdd0);
d2Tnn0 = d2Tdth2 (qq(7:9),Tnn0,dTnn0);

Ty0g = TFromTheta (PP(1:3));
Td0g = TFromTheta (PP(4:6));
Tn0d = TFromTheta (PP(7:9));

Tny0 = (Ty0g.')*Td0g*Tdd0*Tn0d*Tnn0;
Td0r = ((Ty0g*Tyy0*Try).')*Td0g;
Tn0r = Td0r*Tdd0*Tn0d;

Tnr  = Tn0r*Tnn0;
rx   = -Tnr(1,2);
ry   = -Tnr(2,2);
r2   = rx^2 + ry^2;
r22  = r2^2;
azi  = atan2c (ry,rx);

dpdx    = -ry/r2;
dpdy    =  rx/r2;
d2pdx2  =  2*rx*ry/r22;
d2pdy2  = -d2pdx2;
d2pdxdy =  (ry^2 - rx^2)/r22;

drdq = zeros(2,9);
for jj = 1:3

   ic3 = 3*(jj-1);

   dT = ((dTyy0(:,ic3+[1:3])*Try).')*Tny0;
   drdq(:,jj)   = -dT(1:2,2);

   dT = Td0r*dTdd0(:,ic3+[1:3])*Tn0d*Tnn0;
   drdq(:,jj+3) = -dT(1:2,2);

   dT = Tn0r*dTnn0(:,ic3+[1:3]);
   drdq(:,jj+6) = -dT(1:2,2);

end

dpdq = dpdx*drdq(1,:) + dpdy*drdq(2,:);
Dy(1,dofs) = dpdq;

W = dpdq*dqq;

d2rdq2 = zeros(2,81);

for ii = 1:3

   ic9 = 9*(ii-1);
   ic3 = 3*(ii-1);

   for jj = 1:3

      jc3 = 3*(jj-1);

      % q1:3, q1:3
      d2T = ((d2Tyy0(:,ic9+jc3+[1:3])*Try).')*Tny0;
      d2rdq2(:,ic9+jj) = -d2T(1:2,2);

      % q1:3, q4:6
      d2T = ((Ty0g*dTyy0(:,ic3+[1:3])*Try).') ...
          * Td0g*dTdd0(:,jc3+[1:3])*Tn0d*Tnn0;
      d2rdq2(:,ic9+3+jj) = -d2T(1:2,2);

      % q1:3, q7:9
      d2T = ((Ty0g*dTyy0(:,ic3+[1:3])*Try).') ...
          * Td0g*Tdd0*Tn0d*dTnn0(:,jc3+[1:3]);
      d2rdq2(:,ic9+6+jj) = -d2T(1:2,2);

      % q4:6, q1:3
      d2T = ((Ty0g*dTyy0(:,jc3+[1:3])*Try).') ...
          * Td0g*dTdd0(:,ic3+[1:3])*Tn0d*Tnn0;
      d2rdq2(:,27+ic9+jj) = -d2T(1:2,2);

      % q4:6, q4:6
      d2T = Td0r*d2Tdd0(:,ic9+jc3+[1:3])*Tn0d*Tnn0;
      d2rdq2(:,27+ic9+3+jj) = -d2T(1:2,2);

      % q4:6, q7:9
      d2T = Td0r*dTdd0(:,ic3+[1:3])*Tn0d*dTnn0(:,jc3+[1:3]);
      d2rdq2(:,27+ic9+6+jj) = -d2T(1:2,2);

      % q7:9, q1:3
      d2T = ((Ty0g*dTyy0(:,jc3+[1:3])*Try).') ...
          * Td0g*Tdd0*Tn0d*dTnn0(:,ic3+[1:3]);
      d2rdq2(:,54+ic9+jj) = -d2T(1:2,2);

      % q7:9, q4:6
      d2T = Td0r*dTdd0(:,jc3+[1:3])*Tn0d*dTnn0(:,ic3+[1:3]);
      d2rdq2(:,54+ic9+3+jj) = -d2T(1:2,2);

      % q7:9, q7:9
      d2T = Td0r*Tdd0*Tn0d*d2Tnn0(:,ic9+jc3+[1:3]);
      d2rdq2(:,54+ic9+6+jj) = -d2T(1:2,2);

   end

end

rq = drdq*dqq;
dWdqd = dpdq;
dWdq  = d2pdx2*rq(1)*drdq(1,:)   ...
      + d2pdy2*rq(2)*drdq(2,:)   ...
      + d2pdxdy*rq(1)*drdq(2,:)  ...
      + d2pdxdy*rq(2)*drdq(1,:);
for pp = 1:9
   ic9 = 9*(pp-1);
   dWdq(pp) = dWdq(pp) ...
            + dpdx*(d2rdq2(1,ic9+[1:9])*dqq) ...
            + dpdy*(d2rdq2(2,ic9+[1:9])*dqq);
end

Dy(2,dofs) = dWdq;
Dy(2,Ndj+dofs) = dWdqd;
