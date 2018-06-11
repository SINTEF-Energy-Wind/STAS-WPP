%function Lambda = constraints (qB,PB,Tbh,idofs,idofm)
function [Lambda,L,C,ret,slv] = constraints (qB,PB,Tbh,idofs,idofm)
%
% Build the constraint matrix L, then partition and compute Lambda, which
% is used to reduce the equations of motion.  After introducing redundant
% (and covenient) degrees-of-freedom at each joint, there are 6x6 + 3 = 39 
% slave degrees-of-freedom, 6 each at the tower-nacelle, rear bearing,
% and three pitch bearings, as well as 3 positional DOFs at the front 
% bearing.
%
% Version:        Changes:
% --------        -------------
% 21.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 21.10.2017      L = dC/dq verified by complex step.  Equations for C
%                 verified by independent re-coding in a test program,
%                 and comparison with C output from this function, for
%                 a case of an operating wind turbine.
%
% Inputs:
% -------
% qB              : Nodal DOFs augmented with yaw, azimuth, front bearing
%                   axial offset, and three blade pitch DOFs.
% PB              : Nodal offsets and initial body orientation.
% Tbh             : Blade-to-hub coordinate transform.
% idofs,idofm     : Joint DOF references.
%
% Outputs:
% --------
% Lambda          : Partitioned and inverted constraint matrix.
% L               : Matrix of constraint equation variations, L(q) dq = 0.
% C               : Vector of constraint equation values.
% ret,slv         : Retained and slave DOFs.

Ndof = size(qB,1) - 6;  % Let Ndof = 6*Nnod in the unconstrained structure.

Nslv = 39;
L = spalloc (Nslv,Ndof+6,12*39);
C = zeros(Nslv,1);

ijnt = Ndof + [1:6];

yaw0 = PB(ijnt(1));
azi0 = PB(ijnt(2));
fbx0 = PB(ijnt(3));
bt10 = PB(ijnt(4));
bt20 = PB(ijnt(5));
bt30 = PB(ijnt(6));

cy0  = cos(yaw0);
sy0  = sin(yaw0);
ca0  = cos(azi0);
sa0  = sin(azi0);
cb10 = cos(bt10);
sb10 = sin(bt10);
cb20 = cos(bt20);
sb20 = sin(bt20);
cb30 = cos(bt30);
sb30 = sin(bt30);

yaw  = qB(ijnt(1));
azi  = qB(ijnt(2));
fbx  = qB(ijnt(3));
bt1  = qB(ijnt(4));
bt2  = qB(ijnt(5));
bt3  = qB(ijnt(6));

cy   = cos(yaw + yaw0);
sy   = sin(yaw + yaw0);
ca   = cos(azi + azi0);
sa   = sin(azi + azi0);
cb1  = cos(bt1 + bt10);
sb1  = sin(bt1 + bt10);
cb2  = cos(bt2 + bt20);
sb2  = sin(bt2 + bt20);
cb3  = cos(bt3 + bt30);
sb3  = sin(bt3 + bt30);

% Foundation-tower joint.  This is rigid, so all the tower reference
% (slave) DOFs are constrained.
TF0g         = TFromTheta (PB(idofs(1)+[4:6]));
[TFF0,dTFF0] = dTdth      (qB(idofs(1)+[4:6]));
Tm0F         = TFromTheta (PB(idofm(2)+[4:6]));
[Tmm0,dTmm0] = dTdth      (qB(idofm(2)+[4:6]));
TTm          = [0 0 1;1 0 0;0 1 0];
TT0g         = TFromTheta (PB(idofs(2)+[4:6]));
[TTT0,dTTT0] = dTdth      (qB(idofs(2)+[4:6]));
TTg_T        = TT0g*TTT0;
TTg_F        = TF0g*TFF0*Tm0F*Tmm0*TTm;  % Necessary for complex step.
TgT          = TTg_T.';
TgF          = Tm0F*Tmm0*TTm*TgT;
Tm0g         = TF0g*TFF0*Tm0F;
Tgm          = TTm*TgT;

ir = 0;
C(ir+[1:3]) = TF0g*TFF0*(PB(idofm(2)+[1:3]) + qB(idofm(2)+[1:3])) ...
            + PB(idofs(1)+[1:3]) + qB(idofs(1)+[1:3])             ...
            - PB(idofs(2)+[1:3]) - qB(idofs(2)+[1:3]);
L(ir+[1:3],idofs(1)+[1:3]) = eye(3);
L(ir+[1:3],idofm(2)+[1:3]) = TF0g*TFF0;
for idof = 1:3
   ic3 = 3*(idof-1);
   L(ir+[1:3],idofs(1)+idof+3) = TF0g*dTFF0(:,ic3+[1:3]) ...
                               * (PB(idofm(2)+[1:3]) + qB(idofm(2)+[1:3]));
end
L(ir+[1:3],idofs(2)+[1:3]) = -eye(3);

ir = ir + 3;
mat = TTg_F*TgT;                % <== order here is very important.
C(ir+[1:3]) = spinToVecU (mat);
for idof = 1:3
   ic3 = 3*(idof-1);

   S = TF0g*dTFF0(:,ic3+[1:3])*TgF;
   v = spinToVecU (S);
   L(ir+[1:3],idofs(1)+idof+3) = v;

   S = Tm0g*dTmm0(:,ic3+[1:3])*Tgm;
   v = spinToVecU (S);
   L(ir+[1:3],idofm(2)+idof+3) = v;

   S = TTg_F*((TT0g*dTTT0(:,ic3+[1:3])).');
   v = spinToVecU (S);
   L(ir+[1:3],idofs(2)+idof+3) = v;

end

% Tower-nacelle joint.
Tm0T         = TFromTheta (PB(idofm(3)+[4:6]));
[Tmm0,dTmm0] = dTdth      (qB(idofm(3)+[4:6]));
Tym          = [0 0 1;1 0 0;0 1 0]*[cy -sy 0;sy cy 0;0 0 1];
dTym         = [0 0 1;1 0 0;0 1 0]*[-sy -cy 0;cy -sy 0;0 0 0];
Ty0g         = TFromTheta (PB(idofs(3)+[4:6]));
[Tyy0,dTyy0] = dTdth      (qB(idofs(3)+[4:6]));
Tyg_Y        = Ty0g*Tyy0;
Tyg_T        = TT0g*TTT0*Tm0T*Tmm0*Tym;  % Necessary for complex step.
Tgy          = Tyg_Y.';
Tm0g         = TT0g*TTT0*Tm0T;
Tgm          = Tym*Tgy;
Tmg          = TT0g*TTT0*Tm0T*Tmm0;      % Necessary for complex step.

ir = ir + 3;
C(ir+[1:3]) = TT0g*TTT0*(PB(idofm(3)+[1:3]) + qB(idofm(3)+[1:3])) ...
            + PB(idofs(2)+[1:3]) + qB(idofs(2)+[1:3])             ...
            - PB(idofs(3)+[1:3]) - qB(idofs(3)+[1:3]);

L(ir+[1:3],idofs(2)+[1:3]) = eye(3);
L(ir+[1:3],idofm(3)+[1:3]) = TT0g*TTT0;
for idof = 1:3
   ic3 = 3*(idof-1);
   L(ir+[1:3],idofs(2)+idof+3) = TT0g*dTTT0(:,ic3+[1:3]) ...
                               * (PB(idofm(3)+[1:3]) + qB(idofm(3)+[1:3]));
end
L(ir+[1:3],idofs(3)+[1:3]) = -eye(3);

ir = ir + 3;
mat = Tyg_T*Tgy;
C(ir+[1:3]) = spinToVecU (mat);

for idof = 1:3
   ic3 = 3*(idof-1);

   S = TT0g*dTTT0(:,ic3+[1:3])*Tm0T*Tmm0*Tgm;
   v = spinToVecU (S);
   L(ir+[1:3],idofs(2)+idof+3) = v;

   S = Tm0g*dTmm0(:,ic3+[1:3])*Tgm;
   v = spinToVecU (S);
   L(ir+[1:3],idofm(3)+idof+3) = v;

   S = Tyg_T*((Ty0g*dTyy0(:,ic3+[1:3])).');
   v = spinToVecU (S);
   L(ir+[1:3],idofs(3)+idof+3) = v;

end
S = Tmg*dTym*Tgy;
v = spinToVecU (S);
L(ir+[1:3],ijnt(1)) = v;

%{
ijnt(1)
Tmg
dTym
Tgy
Tmg*dTym*Tgy
L(ir+[1:3],ijnt(1))
'---------------------'
%}

% Rear bearing.
Tm0y         = TFromTheta (PB(idofm(4)+[4:6]));
[Tmm0,dTmm0] = dTdth      (qB(idofm(4)+[4:6]));
Tdm          = [0 0 -1;-1 0 0;0 1 0]*[ca -sa 0;sa ca 0;0 0 1];
dTdm         = [0 0 -1;-1 0 0;0 1 0]*[-sa -ca 0;ca -sa 0;0 0 0];
Td0g         = TFromTheta (PB(idofs(4)+[4:6]));
[Tdd0,dTdd0] = dTdth      (qB(idofs(4)+[4:6]));
Tdg_D        = Td0g*Tdd0;
Tdg_Y        = Ty0g*Tyy0*Tm0y*Tmm0*Tdm;  % Necessary for complex step.
Tgd          = Tdg_D.';
Tm0g         = Ty0g*Tyy0*Tm0y;
Tgm          = Tdm*Tgd;
Tmg          = Ty0g*Tyy0*Tm0y*Tmm0;      % Necessary for complex step.

ir = ir + 3;
C(ir+[1:3]) = Ty0g*Tyy0*(PB(idofm(4)+[1:3]) + qB(idofm(4)+[1:3])) ...
            + PB(idofs(3)+[1:3]) + qB(idofs(3)+[1:3])             ...
            - PB(idofs(4)+[1:3]) - qB(idofs(4)+[1:3]);

L(ir+[1:3],idofs(3)+[1:3]) = eye(3);
L(ir+[1:3],idofm(4)+[1:3]) = Ty0g*Tyy0;
for idof = 1:3
   ic3 = 3*(idof-1);
   L(ir+[1:3],idofs(3)+idof+3) = Ty0g*dTyy0(:,ic3+[1:3]) ...
                               * (PB(idofm(4)+[1:3]) + qB(idofm(4)+[1:3]));
end
L(ir+[1:3],idofs(4)+[1:3]) = -eye(3);

ir = ir + 3;
mat = Tdg_Y*Tgd;
C(ir+[1:3]) = spinToVecU (mat);
for idof = 1:3
   ic3 = 3*(idof-1);

   S = Ty0g*dTyy0(:,ic3+[1:3])*Tm0y*Tmm0*Tgm;
   v = spinToVecU (S);
   L(ir+[1:3],idofs(3)+idof+3) = v;

   S = Tm0g*dTmm0(:,ic3+[1:3])*Tgm;
   v = spinToVecU (S);
   L(ir+[1:3],idofm(4)+idof+3) = v;

   S = Tdg_Y*((Td0g*dTdd0(:,ic3+[1:3])).');
   v = spinToVecU (S);
   L(ir+[1:3],idofs(4)+idof+3) = v;

end
S = Tmg*dTdm*Tgd;
v = spinToVecU (S);
L(ir+[1:3],ijnt(2)) = v;

% Front bearing.
Tm0y         = TFromTheta (PB(idofm(5)+[4:6]));
[Tmm0,dTmm0] = dTdth      (qB(idofm(5)+[4:6]));

ir = ir + 3;
C(ir+[1:3]) = Ty0g*Tyy0*(PB(idofm(5)+[1:3]) + qB(idofm(5)+[1:3])  ...
            +            Tm0y*Tmm0*[fbx+fbx0;0;0])                ...
            - Td0g*Tdd0*(PB(idofs(5)+[1:3]) + qB(idofs(5)+[1:3])) ...
            + PB(idofs(3)+[1:3]) + qB(idofs(3)+[1:3])             ...
            - PB(idofs(4)+[1:3]) - qB(idofs(4)+[1:3]);
L(ir+[1:3],idofs(3)+[1:3]) = eye(3);
L(ir+[1:3],idofm(5)+[1:3]) = Tyg_Y;
L(ir+[1:3],idofs(4)+[1:3]) = -eye(3);
L(ir+[1:3],idofs(5)+[1:3]) = -Tdg_D;
L(ir+[1:3],ijnt(3)) = Tyg_Y*Tm0y*Tmm0*[1;0;0];
for idof = 1:3
   ic3 = 3*(idof-1);
   L(ir+[1:3],idofs(3)+idof+3) = Ty0g*dTyy0(:,ic3+[1:3])                  ...
                               * (PB(idofm(5)+[1:3]) + qB(idofm(5)+[1:3]) ...
                               +  Tm0y*Tmm0*[fbx+fbx0;0;0]);
   L(ir+[1:3],idofm(5)+idof+3) = Tyg_Y*Tm0y*dTmm0(:,ic3+[1:3])*[fbx+fbx0;0;0];
   L(ir+[1:3],idofs(4)+idof+3) = -Td0g*dTdd0(:,ic3+[1:3]) ...
                               * (PB(idofs(5)+[1:3]) + qB(idofs(5)+[1:3]));
end

% Blade 1.
Tm0d         = TFromTheta (PB(idofm(6)+[4:6]));
[Tmm0,dTmm0] = dTdth      (qB(idofm(6)+[4:6]));
Tpm          = Tbh*[1 0 0;0 cb1 sb1;0 -sb1 cb1];
dTpm         = Tbh*[0 0 0;0 -sb1 cb1;0 -cb1 -sb1];
Tp0g         = TFromTheta (PB(idofs(6)+[4:6]));
[Tpp0,dTpp0] = dTdth      (qB(idofs(6)+[4:6]));
Tpg_P        = Tp0g*Tpp0;
Tpg_D        = Td0g*Tdd0*Tm0d*Tmm0*Tpm;  % Necessary for complex step.
Tgp          = Tpg_P.';
Tm0g         = Td0g*Tdd0*Tm0d;
Tgm          = Tpm*Tgp;
Tmg          = Td0g*Tdd0*Tm0d*Tmm0;      % Necessary for complex step.

ir = ir + 3;
C(ir+[1:3]) = Td0g*Tdd0*(PB(idofm(6)+[1:3]) + qB(idofm(6)+[1:3])) ...
            + PB(idofs(4)+[1:3]) + qB(idofs(4)+[1:3])             ...
            - PB(idofs(6)+[1:3]) - qB(idofs(6)+[1:3]);
L(ir+[1:3],idofs(4)+[1:3]) = eye(3);
L(ir+[1:3],idofm(6)+[1:3]) = Td0g*Tdd0;
for idof = 1:3
   ic3 = 3*(idof-1);
   L(ir+[1:3],idofs(4)+idof+3) = Td0g*dTdd0(:,ic3+[1:3]) ...
                               * (PB(idofm(6)+[1:3]) + qB(idofm(6)+[1:3]));
end
L(ir+[1:3],idofs(6)+[1:3]) = -eye(3);

ir = ir + 3;
mat = Tpg_D*Tgp;
C(ir+[1:3]) = spinToVecU (mat);
for idof = 1:3
   ic3 = 3*(idof-1);

   S = Td0g*dTdd0(:,ic3+[1:3])*Tm0d*Tmm0*Tgm;
   v = spinToVecU (S);
   L(ir+[1:3],idofs(4)+idof+3) = v;

   S = Tm0g*dTmm0(:,ic3+[1:3])*Tgm;
   v = spinToVecU (S);
   L(ir+[1:3],idofm(6)+idof+3) = v;

   S = Tpg_D*((Tp0g*dTpp0(:,ic3+[1:3])).');
   v = spinToVecU (S);
   L(ir+[1:3],idofs(6)+idof+3) = v;

end
S = Tmg*dTpm*Tgp;
v = spinToVecU (S);
L(ir+[1:3],ijnt(4)) = v;

% Blade 2.
Tm0d         = TFromTheta (PB(idofm(7)+[4:6]));
[Tmm0,dTmm0] = dTdth      (qB(idofm(7)+[4:6]));
Tpm          = Tbh*[1 0 0;0 cb2 sb2;0 -sb2 cb2];
dTpm         = Tbh*[0 0 0;0 -sb2 cb2;0 -cb2 -sb2];
Tp0g         = TFromTheta (PB(idofs(7)+[4:6]));
[Tpp0,dTpp0] = dTdth      (qB(idofs(7)+[4:6]));
Tpg_P        = Tp0g*Tpp0;
Tpg_D        = Td0g*Tdd0*Tm0d*Tmm0*Tpm;  % Necessary for complex step.
Tgp          = Tpg_P.';
Tm0g         = Td0g*Tdd0*Tm0d;
Tgm          = Tpm*Tgp;
Tmg          = Td0g*Tdd0*Tm0d*Tmm0;      % Necessary for complex step.

ir = ir + 3;
C(ir+[1:3]) = Td0g*Tdd0*(PB(idofm(7)+[1:3]) + qB(idofm(7)+[1:3])) ...
            + PB(idofs(4)+[1:3]) + qB(idofs(4)+[1:3])             ...
            - PB(idofs(7)+[1:3]) - qB(idofs(7)+[1:3]);
L(ir+[1:3],idofs(4)+[1:3]) = eye(3);
L(ir+[1:3],idofm(7)+[1:3]) = Td0g*Tdd0;
for idof = 1:3
   ic3 = 3*(idof-1);
   L(ir+[1:3],idofs(4)+idof+3) = Td0g*dTdd0(:,ic3+[1:3]) ...
                               * (PB(idofm(7)+[1:3]) + qB(idofm(7)+[1:3]));
end
L(ir+[1:3],idofs(7)+[1:3]) = -eye(3);

ir = ir + 3;
mat = Tpg_D*Tgp;
C(ir+[1:3]) = spinToVecU (mat);
for idof = 1:3
   ic3 = 3*(idof-1);

   S = Td0g*dTdd0(:,ic3+[1:3])*Tm0d*Tmm0*Tgm;
   v = spinToVecU (S);
   L(ir+[1:3],idofs(4)+idof+3) = v;

   S = Tm0g*dTmm0(:,ic3+[1:3])*Tgm;
   v = spinToVecU (S);
   L(ir+[1:3],idofm(7)+idof+3) = v;

   S = Tpg_D*((Tp0g*dTpp0(:,ic3+[1:3])).');
   v = spinToVecU (S);
   L(ir+[1:3],idofs(7)+idof+3) = v;

end
S = Tmg*dTpm*Tgp;
v = spinToVecU (S);
L(ir+[1:3],ijnt(5)) = v;

% Blade 3.
Tm0d         = TFromTheta (PB(idofm(8)+[4:6]));
[Tmm0,dTmm0] = dTdth      (qB(idofm(8)+[4:6]));
Tpm          = Tbh*[1 0 0;0 cb3 sb3;0 -sb3 cb3];
dTpm         = Tbh*[0 0 0;0 -sb3 cb3;0 -cb3 -sb3];
Tp0g         = TFromTheta (PB(idofs(8)+[4:6]));
[Tpp0,dTpp0] = dTdth      (qB(idofs(8)+[4:6]));
Tpg_P        = Tp0g*Tpp0;
Tpg_D        = Td0g*Tdd0*Tm0d*Tmm0*Tpm;  % Necessary for complex step.
Tgp          = Tpg_P.';
Tm0g         = Td0g*Tdd0*Tm0d;
Tgm          = Tpm*Tgp;
Tmg          = Td0g*Tdd0*Tm0d*Tmm0;      % Necessary for complex step.

ir = ir + 3;
C(ir+[1:3]) = Td0g*Tdd0*(PB(idofm(8)+[1:3]) + qB(idofm(8)+[1:3])) ...
            + PB(idofs(4)+[1:3]) + qB(idofs(4)+[1:3])             ...
            - PB(idofs(8)+[1:3]) - qB(idofs(8)+[1:3]);
L(ir+[1:3],idofs(4)+[1:3]) = eye(3);
L(ir+[1:3],idofm(8)+[1:3]) = Td0g*Tdd0;
for idof = 1:3
   ic3 = 3*(idof-1);
   L(ir+[1:3],idofs(4)+idof+3) = Td0g*dTdd0(:,ic3+[1:3]) ...
                               * (PB(idofm(8)+[1:3]) + qB(idofm(8)+[1:3]));
end
L(ir+[1:3],idofs(8)+[1:3]) = -eye(3);

ir = ir + 3;
mat = Tpg_D*Tgp;
C(ir+[1:3]) = spinToVecU (mat);
for idof = 1:3
   ic3 = 3*(idof-1);

   S = Td0g*dTdd0(:,ic3+[1:3])*Tm0d*Tmm0*Tgm;
   v = spinToVecU (S);
   L(ir+[1:3],idofs(4)+idof+3) = v;

   S = Tm0g*dTmm0(:,ic3+[1:3])*Tgm;
   v = spinToVecU (S);
   L(ir+[1:3],idofm(8)+idof+3) = v;

   S = Tpg_D*((Tp0g*dTpp0(:,ic3+[1:3])).');
   v = spinToVecU (S);
   L(ir+[1:3],idofs(8)+idof+3) = v;

end
S = Tmg*dTpm*Tgp;
v = spinToVecU (S);
L(ir+[1:3],ijnt(6)) = v;

slv = slaveDOFs(idofs);
[Lp,rr,ret] = partitionMatrix(L,[],slv);
Nret = size(ret,1);
Nslv = size(slv,1);

Lambda = sparse([eye(Nret);-(Lp(:,Nret+[1:Nslv])\Lp(:,[1:Nret]))]);

%cond(Lp(:,Nret+[1:Nslv]))
%Lp(:,Nret+[1:Nslv])
%Lp(:,[1:Nret])

