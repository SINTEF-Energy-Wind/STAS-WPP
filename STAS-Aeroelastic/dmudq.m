function dmu = dmudq (qn1,qn2,Pn1,Pn2,TsB,dTsB)
%
% mu includes delta_x, the elastic stretch of the element, as well
% as zeta, the elastic portion of the nodal rotation, associated
% with a given element.  The stretch is
%
% (r_k+1^s)x - (r_k^s)_x - L;  L = r0_k+1^s0 - r0_k^s0.
%
% The elastic rotation is
%
%                --                     --
%              1 | (T_n^s)32 - (T_n^s)23 |
% zeta_n\s^s = - | (T_n^s)13 - (T_n^s)31 |
%              2 | (T_n^s)21 - (T_n^s)12 |
%                --                     --
%
% T_n^s = T_B^s*T_s0^B*T_n^n0
%
% dT_n^s   dT_B^s                              dT_n^n0
% ------ = ------*T_s0^B*T_n^n0 + T_B^s*T_s0^B*-------
%   dq       dq                                  dq
%
% Version:        Changes:
% --------        -------------
% 11.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 11.10.2017      Verified by complex step and finite difference.
%
% Inputs:
% -------
% qn1,qn2         : 6*Nel vectors for node 1 and 2 DOFs.
% Pn1,Pn2         : 6*Nel vectors describing undeformed nodal position
%                   and orientation of nodes.
% TsB             : 3-by-iel set of transforms from section to body
%                   coordinates.
% dTsB            : A 3-by-3*12*Nel matrix, derivatives wrt xn1, xn2.
%
% Outputs:
% --------
% dmu             : dmu/dq.
%

Nel = size(qn1,1)/6;

dmu = zeros (12,12*Nel);

for iel = 1:Nel

   ir3  = 3*(iel-1);
   ir6  = 6*(iel-1);
   ir12 = 12*(iel-1);
   ir36 = 36*(iel-1);

   Ts0B           = TFromTheta (Pn2(ir6+[4:6]));

   Tnn0           = zeros(3,6);
   dTnn0          = zeros(3,18);

   [T,dT]         = dTdth (qn1(ir6+[4:6]));
   Tnn0(:,1:3)    = T;
   dTnn0(:,1:9)   = dT;
   [T,dT]         = dTdth (qn2(ir6+[4:6]));
   Tnn0(:,4:6)    = T;
   dTnn0(:,10:18) = dT;

   xn1 = zeros(6,1);
   xn2 = zeros(6,1);
   xn1(1:3) = Pn1(ir6+[1:3]) + qn1(ir6+[1:3]);
   xn1(4:6) = qn1(ir6+[4:6]);
   xn2(1:3) = Pn2(ir6+[1:3]) + qn2(ir6+[1:3]);
   xn2(4:6) = qn2(ir6+[4:6]);

   dx = xn2(1:3) - xn1(1:3);

   Ts0jB         = TFromTheta (Pn2(ir6+[4:6]));
   Ts0j1B        = TFromTheta (Pn1(ir6+[4:6]));
   Tn0ks0j       = (Ts0jB.')*Ts0j1B;
   Ts0jn0k       = Tn0ks0j.';

   for ix = 1:3

      ix3  = 3*(ix-1);

      ee = zeros(3,1);
      ee(ix) = 1;

      dmu(1,ir12+ix) = -[0.5 0 0]*(-(TsB(:,ir3+[1:3]).')*ee ...
                     +             (dTsB(:,ir36+ix3+[1:3]).')*dx);

      dmu(7,ir12+ix) = -dmu(1,ir12+ix);

      dTns = (dTsB(:,ir36+ix3+[1:3]).')*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k;
      dmu(4:6,ir12+ix) = spinToVec(dTns);

      dTns = (dTsB(:,ir36+ix3+[1:3]).')*Ts0B*Tnn0(:,4:6);
      dmu(10:12,ir12+ix) = spinToVec(dTns);

%{
'-------------------------------------------------'
ix
dTsB(:,ir36+ix3+[1:3])
Ts0B
Tnn0(:,1:3)
dTsB(:,ir36+ix3+[1:3])*Ts0B*Tnn0(:,1:3)
%}

   end

   for ix = 4:6

      ix3  = 3*(ix-1);
      in3  = 3*(ix-4);

      dmu(1,ir12+ix) = -[0.5 0 0]*(dTsB(:,ir36+ix3+[1:3]).')*dx;
      dmu(7,ir12+ix) = -dmu(1,ir12+ix);

      dTns = (dTsB(:,ir36+ix3+[1:3]).')*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k ...
           + (TsB(:,ir3+[1:3]).')*Ts0B*Tn0ks0j*dTnn0(:,in3+[1:3])*Ts0jn0k;
      dmu(4:6,ir12+ix) = spinToVec(dTns);

%{
if (imag(qn1(5)) ~= 0)
'dmudq.m'
ix
del = sqrt(eps);
d2Tnsc = imag(dTns)/del;
'dmudq term 1'
dTsB(:,ir36+ix3+[1:3]).'
Ts0B
imag(Tnn0(:,1:3))
(dTsB(:,ir36+ix3+[1:3]).')*Ts0B*Tnn0(:,1:3)
'dmudq term 2'
TsB(:,ir3+[1:3]).'
Ts0B
imag(dTnn0(:,in3+[1:3]))
(TsB(:,ir3+[1:3]).')*Ts0B*dTnn0(:,in3+[1:3])
d2Tnsc
end
%}

      dTns = (dTsB(:,ir36+ix3+[1:3]).')*Ts0B*Tnn0(:,4:6);
      dmu(10:12,ir12+ix) = spinToVec(dTns);

%{
'-------------------------------------------------'
ix
full(dTsB(:,ir36+ix3+[1:3]))
full(TsB)
Ts0B
Tnn0
dTnn0(:,in3+[1:3])

TsB(ir3+[1:3])
Ts0B
dTnn0(:,in3+[1:3])
TsB(ir3+[1:3])*Ts0B*dTnn0(:,in3+[1:3])

dTsB(:,ir36+ix3+[1:3])*Ts0B*Tnn0(:,1:3) ...
+ TsB(ir3+[1:3])*Ts0B*dTnn0(:,in3+[1:3])
dTsB(:,ir36+ix3+[1:3])*Ts0B*Tnn0(:,4:6)
%}


   end

   for ix = 7:9

      ix3  = 3*(ix-1);

      ee = zeros(3,1);
      ee(ix-6) = 1;

      dmu(1,ir12+ix) = -[0.5 0 0]*((TsB(:,ir3+[1:3]).')*ee ...
                     +             (dTsB(:,ir36+ix3+[1:3]).')*dx);

      dmu(7,ir12+ix) = -dmu(1,ir12+ix);

      dTns = (dTsB(:,ir36+ix3+[1:3]).')*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k;
      dmu(4:6,ir12+ix) = spinToVec(dTns);

      dTns = (dTsB(:,ir36+ix3+[1:3]).')*Ts0B*Tnn0(:,4:6);
      dmu(10:12,ir12+ix) = spinToVec(dTns);

   end

   for ix = 10:12

      ix3  = 3*(ix-1);
      in3  = 3*(ix-10);

      dmu(1,ir12+ix) = -[0.5 0 0]*(dTsB(:,ir36+ix3+[1:3]).')*dx;
      dmu(7,ir12+ix) = -dmu(1,ir12+ix);

      dTns = (dTsB(:,ir36+ix3+[1:3]).')*Ts0B*Tn0ks0j*Tnn0(:,1:3)*Ts0jn0k;
      dmu(4:6,ir12+ix) = spinToVec(dTns);

      dTns = (dTsB(:,ir36+ix3+[1:3]).')*Ts0B*Tnn0(:,4:6) ...
           + (TsB(:,ir3+[1:3]).')*Ts0B*dTnn0(:,in3+9+[1:3]);
      dmu(10:12,ir12+ix) = spinToVec(dTns);

   end

end
