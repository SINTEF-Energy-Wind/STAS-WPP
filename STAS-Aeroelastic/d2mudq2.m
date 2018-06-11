function d2mu = d2mudq2 (qn1,qn2,Pn1,Pn2,TsB,dTsB,d2TsB)
%
% Take the second derivative of mu, which consists of the elastic stretch
% and elastic nodal rotations.
%
% Version:        Changes:
% --------        -------------
% 13.11.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 13.11.2017      Verified by complex step derivatives using dmudq.m.
%
% Inputs:
% -------
% qn1,qn2         : 6*Nel vectors for node 1 and 2 DOFs.
% Pn1,Pn2         : 6*Nel vectors describing undeformed nodal position
%                   and orientation of nodes.
% TsB             : 3-by-3 transform from section to body coordinates.
% dTsB            : A 3-by-3*12 matrix, derivatives wrt xn1, xn2.
% d2TsB           : Second derivatives.  3-by-3*12*12.
%
% Outputs:
% --------
% d2mu            : d^2mu/dq^2.

d2mu = zeros(12,12*12);

aB   = Pn2(1:3) + qn2(1:3) - Pn1(1:3) - qn1(1:3);

TBs  = TsB.';
Ts0B = TFromTheta (Pn2(4:6));
TBs0 = Ts0B.';

Tnn0            = zeros(3,6);
dTnn0           = zeros(3,3*6);
d2Tnn0          = zeros(3,3*3*6);
[T,dT]          = dTdth (qn1(4:6));
d2T             = d2Tdth2 (qn1(4:6),T,dT);
Tnn0(:,1:3)     = T;
dTnn0(:,1:9)    = dT;
d2Tnn0(:,1:27)  = d2T;
[T,dT]          = dTdth (qn2(4:6));
d2T             = d2Tdth2 (qn2(4:6),T,dT);
Tnn0(:,4:6)     = T;
dTnn0(:,10:18)  = dT;
d2Tnn0(:,28:54) = d2T;

Ts0jB         = TFromTheta (Pn2(4:6));
Ts0j1B        = TFromTheta (Pn1(4:6));
Tn0ks0j       = (Ts0jB.')*Ts0j1B;
Ts0jn0k       = Tn0ks0j.';

% d1 , d1.  Verified.
for ii = 1:3

   ic    =  0;              % d1.
   signi = -1;

   ic36 = 36*(ii+ic-1);
   ic12 = 12*(ii+ic-1);
   ic3  =  3*(ii+ic-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  0;           % d1.
      signj = -1;

      jc3  =  3*(jj+jc-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*(signi*dTBsj(1,ii) + signj*dTBsi(1,jj) ...
                         +       d2TBsij(1,:)*aB);

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% d1 , d2
for ii = 1:3

   ic    =  0;              % d1.
   signi = -1;

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  6;           % d2.
      signj =  1;

      jc3   =  3*(jj+jc-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*(signi*dTBsj(1,ii) + signj*dTBsi(1,jj) ...
                         +       d2TBsij(1,:)*aB);

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% d2 , d1
for ii = 1:3

   ic    =  6;              % d2.
   signi =  1;

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  0;           % d1.
      signj = -1;

      jc3   =  3*(jj+jc-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*(signi*dTBsj(1,ii) + signj*dTBsi(1,jj) ...
                         +       d2TBsij(1,:)*aB);

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% d2 , d2
for ii = 1:3

   ic    =  6;              % d2.
   signi =  1;

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  6;           % d2.
      signj =  1;

      jc3   =  3*(jj+jc-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*(signi*dTBsj(1,ii) + signj*dTBsi(1,jj) ...
                         +       d2TBsij(1,:)*aB);

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% d1 , th1
for ii = 1:3

   ic    =  0;              % d1.
   signi = -1;

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  3;           % th1.

      jc3   =  3*(jj+jc-1);
      j3    =  3*(jj-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*(signi*dTBsj(1,ii) + d2TBsij(1,:)*aB);

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k ...
          + dTBsi*Ts0B*Tn0ks0j*dTnn0(:,j3+[1:3])*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% d1 , th2
for ii = 1:3

   ic    =  0;              % d1.
   signi = -1;

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  9;           % th2.

      jc3   =  3*(jj+jc-1);
      j3    =  3*(jj-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*(signi*dTBsj(1,ii) + d2TBsij(1,:)*aB);

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6) + dTBsi*Ts0B*dTnn0(:,j3+9+[1:3]);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% d2 , th1
for ii = 1:3

   ic    =  6;              % d2.
   signi =  1;

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  3;           % th1.

      jc3   =  3*(jj+jc-1);
      j3    =  3*(jj-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*(signi*dTBsj(1,ii) + d2TBsij(1,:)*aB);

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k ...
          + dTBsi*Ts0B*Tn0ks0j*dTnn0(:,j3+[1:3])*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% d2 , th2
for ii = 1:3

   ic    =  6;              % d2.
   signi =  1;

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  9;           % th2.

      jc3   =  3*(jj+jc-1);
      j3    =  3*(jj-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*(signi*dTBsj(1,ii) + d2TBsij(1,:)*aB);

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6) + dTBsi*Ts0B*dTnn0(:,j3+9+[1:3]);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% th1 , d1
for ii = 1:3

   ic    =  3;              % th1.

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);
   i3    =  3*(ii-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  0;           % d1.
      signj = -1;

      jc3   =  3*(jj+jc-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*(signj*dTBsi(1,jj) + d2TBsij(1,:)*aB);

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k ...
          + dTBsj*Ts0B*Tn0ks0j*dTnn0(:,i3+[1:3])*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% th1 , d2
for ii = 1:3

   ic    =  3;              % th1.

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);
   i3    =  3*(ii-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  6;           % d2.
      signj =  1;

      jc3   =  3*(jj+jc-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*(signj*dTBsi(1,jj) + d2TBsij(1,:)*aB);

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k ...
          + dTBsj*Ts0B*Tn0ks0j*dTnn0(:,i3+[1:3])*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% th2 , d1
for ii = 1:3

   ic    =  9;              % th2.

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);
   i3    =  3*(ii-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  0;           % d1.
      signj = -1;

      jc3   =  3*(jj+jc-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*(signj*dTBsi(1,jj) + d2TBsij(1,:)*aB);

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6) + dTBsj*Ts0B*dTnn0(:,i3+9+[1:3]);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% th2 , d2
for ii = 1:3

   ic    =  9;              % th2.

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);
   i3    =  3*(ii-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  6;           % d2.
      signj =  1;

      jc3   =  3*(jj+jc-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*(signj*dTBsi(1,jj) + d2TBsij(1,:)*aB);

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6) + dTBsj*Ts0B*dTnn0(:,i3+9+[1:3]);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% th1 , th1
for ii = 1:3

   ic    =  3;              % th1.

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);
   i9    =  9*(ii-1);
   i3    =  3*(ii-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  3;           % th1.

      jc3   =  3*(jj+jc-1);
      j3    =  3*(jj-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*d2TBsij(1,:)*aB;

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k      ...
          + dTBsj*Ts0B*Tn0ks0j*dTnn0(:,i3+[1:3])*Ts0jn0k  ...
          + dTBsi*Ts0B*Tn0ks0j*dTnn0(:,j3+[1:3])*Ts0jn0k  ...
          + TBs*Ts0B*Tn0ks0j*d2Tnn0(:,i9+j3+[1:3])*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% th1 , th2
for ii = 1:3

   ic    =  3;              % th1.

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);
   i9    =  9*(ii-1);
   i3    =  3*(ii-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  9;           % th2.

      jc3   =  3*(jj+jc-1);
      j3    =  3*(jj-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*d2TBsij(1,:)*aB;

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k ...
          + dTBsj*Ts0B*Tn0ks0j*dTnn0(:,i3+[1:3])*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6) ...
          + dTBsi*Ts0B*dTnn0(:,j3+9+[1:3]);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% th2 , th1
for ii = 1:3

   ic    =  9;              % th2.

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);
   i9    =  9*(ii-1);
   i3    =  3*(ii-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  3;           % th1.

      jc3   =  3*(jj+jc-1);
      j3    =  3*(jj-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*d2TBsij(1,:)*aB;

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k ...
          + dTBsi*Ts0B*Tn0ks0j*dTnn0(:,j3+[1:3])*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6) ...
          + dTBsj*Ts0B*dTnn0(:,i3+9+[1:3]);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

% th2 , th2
for ii = 1:3

   ic    =  9;              % th2.

   ic36  = 36*(ii+ic-1);
   ic12  = 12*(ii+ic-1);
   ic3   =  3*(ii+ic-1);
   i9    =  9*(ii-1);
   i3    =  3*(ii-1);

   dTBsi   =  dTsB(:,ic3+[1:3]).';

   for jj = 1:3

      jc    =  9;           % th2.

      jc3   =  3*(jj+jc-1);
      j3    =  3*(jj-1);

      dTBsj   =  dTsB(:,jc3+[1:3]).';
      d2TBsij = d2TsB(:,ic36+jc3+[1:3]).';

      d2mu(1,ic12+jc+jj) = -0.5*d2TBsij(1,:)*aB;

      d2mu(7,ic12+jc+jj) = -d2mu(1,ic12+jc+jj);

      d2T = d2TBsij*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k;
      d2mu(4:6,ic12+jc+jj) = spinToVec(d2T);

      d2T = d2TBsij*Ts0B*Tnn0(:,4:6)        ...
          + dTBsj*Ts0B*dTnn0(:,i3+9+[1:3])  ...
          + dTBsi*Ts0B*dTnn0(:,j3+9+[1:3])  ...
          + TBs*Ts0B*d2Tnn0(:,i9+j3+27+[1:3]);
      d2mu(10:12,ic12+jc+jj) = spinToVec(d2T);

   end

end

d2mu = sparse(d2mu);
