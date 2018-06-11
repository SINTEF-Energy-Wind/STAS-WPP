function D = angularVelocity (q,inods,inodm,P0)
%
% The nodal rotational DOFs in STAS consist of an axis-angle
% representation, in the form of an exponential map.  That is, the
% three parameters thx,thy,thz define the orientation of an axis
% with norm(th) the angle of rotation, generating a transformation
% matrix Tb^a.  The axis is defined in "a" coordinates.
%
% The time derivatives of the rotation parameters are related to
% the angular velocity vector by w_b/a^a = D (dth/dt), where D(th) 
% is a matrix that is a function of theta.
%
% To find D, write S = (dT/dt)T', or
% (dT/dthx)T'(dthx/dt) + (dT/dthy)T'(dthy/dt) + (dT/dthz)T'(dthz/dt),
% where S is spin(w_b/a^a).  Then extract the appropriate elements
% of each of the (dT/dthj)T's, corresponding to the w's in S, and
% write these in matrix-vector format.  If we say that Fj = (dT/dthj)T',
% then
%     --                      --
%     | -Fx,23  -Fy,23  -Fz,23 |
% D = |  Fx,13   Fy,13   Fz,13 |
%     | -Fx,12  -Fy,12  -Fz,12 |
%     --                      --
%
% In general, T can be written as T_B0^g*T_B^B0*T_n0^B*T_n^n0.
%
% This takes a vector of DOFs q and outputs the D matrix, which for
% each set of three nodal rotations relates W_B\B0^g and Th_B^B0
% (body reference nodes) or w_n\n0^g and th_n^n0.
%
% w_n\g^g = D1 dTh_B^B0/dt + Dk dth_n^n0/dt
%
% Version:        Changes:
% --------        -------------
% 07.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 07.10.2017      
%
% Inputs:
% -------
% q               : DOF vector.  3 positions, 3 rotations for each node.
% inods,inodm     : Body node indices.
% iels            : Body element indices.
% P0              : Initial undeformed node-to-body positions and
%                   rotations.
%
% Outputs:
% --------
% D              : 3-by-3*Nnod matrix.
%

Nnod = size(q,1)/6;
Nnb  = inods(8) - inods(7);

D = zeros(3,3*Nnod);

inod = 1;
ir3 = 3*(inod-1);
ir6 = 6*(inod-1);
[F,TB_B0,dTB_B0,TB0_g] = FRefMatrix (q(ir6+[4:6]),P0(ir6+[4:6]));
D(:,ir3+[1:3]) = Dmatrix (F);
TBg = TB0_g*TB_B0;

%{
'--------------------------------------------------------------'
[inod P0(ir6+4) P0(ir6+5) P0(ir6+6) q(ir6+4) q(ir6+5) q(ir6+6)]
TB_B0
dTB_B0
TB0_g
TBg
F
D(:,ir3+[1:3])
%}

for inod = inods(1)+1:inodm(2)
   ir3 = 3*(inod-1);
   ir6 = 6*(inod-1);
   Ts0B = TFromTheta (P0(ir6+[4:6]));
   [F,Tnn0,dTnn0] = Fmatrix (q(ir6+[4:6]),TBg,Ts0B);
   D(:,ir3+[1:3]) = Dmatrix (F);

%{
[inod P0(ir6+4) P0(ir6+5) P0(ir6+6) q(ir6+4) q(ir6+5) q(ir6+6)]
TBg
Ts0B
Tnn0
dTnn0
F
D(:,ir3+[1:3])
%}

end

inod = inods(2);
ir3 = 3*(inod-1);
ir6 = 6*(inod-1);
[F,TB_B0,dTB_B0,TB0_g] = FRefMatrix (q(ir6+[4:6]),P0(ir6+[4:6]));
D(:,ir3+[1:3]) = Dmatrix (F);
TBg = TB0_g*TB_B0;

%{
'--------------------------------------------------------------'
[inod P0(ir6+4) P0(ir6+5) P0(ir6+6) q(ir6+4) q(ir6+5) q(ir6+6)]
TB_B0
dTB_B0
TB0_g
TBg
F
D(:,ir3+[1:3])
%}

for inod = inods(2)+1:inodm(3)
   ir3 = 3*(inod-1);
   ir6 = 6*(inod-1);
   Ts0B = TFromTheta (P0(ir6+[4:6]));
   [F,Tnn0,dTnn0] = Fmatrix (q(ir6+[4:6]),TBg,Ts0B);
   D(:,ir3+[1:3]) = Dmatrix (F);

%{
[inod P0(ir6+4) P0(ir6+5) P0(ir6+6) q(ir6+4) q(ir6+5) q(ir6+6)]
TBg
Ts0B
Tnn0
dTnn0
F
D(:,ir3+[1:3])
%}

end

inod = inods(3);
ir3 = 3*(inod-1);
ir6 = 6*(inod-1);
[F,TB_B0,dTB_B0,TB0_g] = FRefMatrix (q(ir6+[4:6]),P0(ir6+[4:6]));
D(:,ir3+[1:3]) = Dmatrix (F);
TBg = TB0_g*TB_B0;

%{
'--------------------------------------------------------------'
[inod P0(ir6+4) P0(ir6+5) P0(ir6+6) q(ir6+4) q(ir6+5) q(ir6+6)]
TB_B0
dTB_B0
TB0_g
TBg
F
D(:,ir3+[1:3])
%}

for inod = inods(3)+1:inodm(5)
   ir3 = 3*(inod-1);
   ir6 = 6*(inod-1);
   Ts0B = TFromTheta (P0(ir6+[4:6]));
   [F,Tnn0,dTnn0] = Fmatrix (q(ir6+[4:6]),TBg,Ts0B);
   D(:,ir3+[1:3]) = Dmatrix (F);

%{
[inod P0(ir6+4) P0(ir6+5) P0(ir6+6) q(ir6+4) q(ir6+5) q(ir6+6)]
TBg
Ts0B
Tnn0
dTnn0
F
D(:,ir3+[1:3])
%}

end

inod = inods(4);
ir3 = 3*(inod-1);
ir6 = 6*(inod-1);
[F,TB_B0,dTB_B0,TB0_g] = FRefMatrix (q(ir6+[4:6]),P0(ir6+[4:6]));
D(:,ir3+[1:3]) = Dmatrix (F);
TBg = TB0_g*TB_B0;

%{
'--------------------------------------------------------------'
[inod P0(ir6+4) P0(ir6+5) P0(ir6+6) q(ir6+4) q(ir6+5) q(ir6+6)]
TB_B0
dTB_B0
TB0_g
TBg
F
D(:,ir3+[1:3])
%}

for inod = inods(4)+1:inodm(8)
   ir3 = 3*(inod-1);
   ir6 = 6*(inod-1);
   Ts0B = TFromTheta (P0(ir6+[4:6]));
   [F,Tnn0,dTnn0] = Fmatrix (q(ir6+[4:6]),TBg,Ts0B);
   D(:,ir3+[1:3]) = Dmatrix (F);

%{
[inod P0(ir6+4) P0(ir6+5) P0(ir6+6) q(ir6+4) q(ir6+5) q(ir6+6)]
TBg
Ts0B
Tnn0
dTnn0
F
D(:,ir3+[1:3])
%}

end

inod = inods(6);
ir3 = 3*(inod-1);
ir6 = 6*(inod-1);
[F,TB_B0,dTB_B0,TB0_g] = FRefMatrix (q(ir6+[4:6]),P0(ir6+[4:6]));
D(:,ir3+[1:3]) = Dmatrix (F);
TBg = TB0_g*TB_B0;

%{
'--------------------------------------------------------------'
[inod P0(ir6+4) P0(ir6+5) P0(ir6+6) q(ir6+4) q(ir6+5) q(ir6+6)]
TB_B0
dTB_B0
TB0_g
TBg
F
D(:,ir3+[1:3])
%}

for inod = inods(6)+[1:Nnb-1]
   ir3 = 3*(inod-1);
   ir6 = 6*(inod-1);
   Ts0B = TFromTheta (P0(ir6+[4:6]));
   [F,Tnn0,dTnn0] = Fmatrix (q(ir6+[4:6]),TBg,Ts0B);
   D(:,ir3+[1:3]) = Dmatrix (F);

%{
[inod P0(ir6+4) P0(ir6+5) P0(ir6+6) q(ir6+4) q(ir6+5) q(ir6+6)]
TBg
Ts0B
Tnn0
dTnn0
F
D(:,ir3+[1:3])
%}

end

inod = inods(7);
ir3 = 3*(inod-1);
ir6 = 6*(inod-1);
[F,TB_B0,dTB_B0,TB0_g] = FRefMatrix (q(ir6+[4:6]),P0(ir6+[4:6]));
D(:,ir3+[1:3]) = Dmatrix (F);
TBg = TB0_g*TB_B0;

%{
'--------------------------------------------------------------'
[inod P0(ir6+4) P0(ir6+5) P0(ir6+6) q(ir6+4) q(ir6+5) q(ir6+6)]
TB_B0
dTB_B0
TB0_g
TBg
F
D(:,ir3+[1:3])
%}

for inod = inods(7)+[1:Nnb-1]
   ir3 = 3*(inod-1);
   ir6 = 6*(inod-1);
   Ts0B = TFromTheta (P0(ir6+[4:6]));
   [F,Tnn0,dTnn0] = Fmatrix (q(ir6+[4:6]),TBg,Ts0B);
   D(:,ir3+[1:3]) = Dmatrix (F);

%{
[inod P0(ir6+4) P0(ir6+5) P0(ir6+6) q(ir6+4) q(ir6+5) q(ir6+6)]
TBg
Ts0B
Tnn0
dTnn0
F
D(:,ir3+[1:3])
%}

end

inod = inods(8);
ir3 = 3*(inod-1);
ir6 = 6*(inod-1);
[F,TB_B0,dTB_B0,TB0_g] = FRefMatrix (q(ir6+[4:6]),P0(ir6+[4:6]));
D(:,ir3+[1:3]) = Dmatrix (F);
TBg = TB0_g*TB_B0;

%{
'--------------------------------------------------------------'
[inod P0(ir6+4) P0(ir6+5) P0(ir6+6) q(ir6+4) q(ir6+5) q(ir6+6)]
TB_B0
dTB_B0
TB0_g
TBg
F
D(:,ir3+[1:3])
%}

for inod = inods(8)+[1:Nnb-1]
   ir3 = 3*(inod-1);
   ir6 = 6*(inod-1);
   Ts0B = TFromTheta (P0(ir6+[4:6]));
   [F,Tnn0,dTnn0] = Fmatrix (q(ir6+[4:6]),TBg,Ts0B);
   D(:,ir3+[1:3]) = Dmatrix (F);

%{
[inod P0(ir6+4) P0(ir6+5) P0(ir6+6) q(ir6+4) q(ir6+5) q(ir6+6)]
TBg
Ts0B
Tnn0
dTnn0
F
D(:,ir3+[1:3])
%}

end

