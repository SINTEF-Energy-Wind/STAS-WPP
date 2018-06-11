function mu = getMu (qn1,qn2,Pn1,Pn2,TsB)
%
% mu includes delta_x, the elastic stretch of the element, as well
% as zeta, the elastic portion of the nodal rotation, associated
% with a given element.  The stretch is
%
% (r_k+1^s)x - (r_k^s)_x - L;  L = |p_k+1^s0 - p_k^s0|.
%
% The elastic rotation is
%
%                --                     --
%              1 | (T_e^s)32 - (T_e^s)23 |
% zeta_n\s^s = - | (T_e^s)13 - (T_e^s)31 |
%              2 | (T_e^s)21 - (T_e^s)12 |
%                --                     --
%
% T_e^s = T_B^s*T_s0^B*T_n0^s0*T_n^n0*T_s0^n0
%
% Version:        Changes:
% --------        -------------
% 11.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 11.10.2017      Checked some simple test cases.  This was also used
%                 in the 45-degree bend case.
%
% Inputs:
% -------
% qn1,qn2         : 6*Nel vectors for node 1 and 2 DOFs.
% Pn1,Pn2         : 6*Nel vectors describing undeformed nodal position
%                   and orientation of nodes.
% TsB             : 3-by-iel set of transforms from section to body
%                   coordinates.
%
% Outputs:
% --------
% mu              : [delta 0 0 zeta]^T.
%

%'getMu'

Nel = size(qn1,1)/6;

mu = zeros (12,Nel);

for iel = 1:Nel

   ir3 = 3*(iel-1);
   ir6 = 6*(iel-1);

   xn1 = zeros(6,1);
   xn2 = zeros(6,1);
   xn1(1:3) = Pn1(ir6+[1:3]) + qn1(ir6+[1:3]);
   xn1(4:6) = qn1(ir6+[4:6]);
   xn2(1:3) = Pn2(ir6+[1:3]) + qn2(ir6+[1:3]);
   xn2(4:6) = qn2(ir6+[4:6]);

   % The undeformed element section transform is given by the
   % second node.
   Ts0jB         = TFromTheta (Pn2(ir6+[4:6]));
   Ts0j1B        = TFromTheta (Pn1(ir6+[4:6]));
   Tn0ks0j       = (Ts0jB.')*Ts0j1B;
   Tnn0          = zeros(3,6);
   Tnn0(:,1:3)   = TFromTheta (qn1(ir6+[4:6]));
   Tnn0(:,4:6)   = TFromTheta (qn2(ir6+[4:6]));

   Teks0j        = Tn0ks0j*Tnn0(:,1:3)*(Tn0ks0j.');
   Tek1s0j       = Tnn0(:,4:6);

   Ts0s          = (TsB(:,ir3+[1:3]).')*Ts0jB;

   Tns           = Ts0s*Teks0j;
   mu(4:6,iel)   = spinToVec (Tns);

   Tns           = Ts0s*Tek1s0j;
   mu(10:12,iel) = spinToVec (Tns);


%[Pn1 qn1 xn1 Pn2 qn2 xn2]

%TsB(:,ir3+[1:3])
%Ts0jB
%(TsB(:,ir3+[1:3]).')
%Ts0jB
%Tn0ks0j
%Tnn0(:,1:3)
%(Tn0ks0j.')
%'--------------'
%TsB(:,ir3+[1:3])
%Ts0jB
%Tnn0(:,4:6)


   L = Pn2(ir6+[1:3]) - Pn1(ir6+[1:3]);
   dx = xn2(1:3) - xn1(1:3);

   mu(1,iel)     = -[0.5 0 0]*((TsB(:,ir3+[1:3]).')*dx ...
                 -             Ts0jB.'*L);
   mu(7,iel)     = -mu(1,iel);

end
