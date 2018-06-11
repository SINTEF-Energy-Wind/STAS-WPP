function dL = derivConstraints (qB,PB,Tbh,idofs,idofm)
%
% The implementation of constraints in the equations of motion includes
% a term representing the centripetal constraint forces.  This term has
% the form Fc = Lambda^(-1)*M*Gamma, where Gamma includes a term
% -Ls^(-1)*dqk/dt*dL/dqk*dq/dt.  This function computes dL/dqk for the
% relevant variables.
%
% Version:        Changes:
% --------        -------------
% 28.11.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 28.11.2017      Derivatives of L verified using complex step, with
%                 constraints.m.
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
% dL              : dL/dq.
%

Ndj  = size(qB,1);
Ndof = Ndj - 6;  % Let Ndof = 6*Nnod in the unconstrained structure.

dL = spalloc (39,72*Ndj,5*Ndj);

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

% Column refs in the dL matrix, for each L submatrix.
ThFref  = Ndj* [0:2].';
dmtref  = Ndj* [3:5].';
thmtref = Ndj* [6:8].';
ThTref  = Ndj* [9:11].';
dmyref  = Ndj*[12:14].';
thmyref = Ndj*[15:17].';
ThYref  = Ndj*[18:20].';
dmdref  = Ndj*[21:23].';
thmdref = Ndj*[24:26].';
ThDref  = Ndj*[27:29].';
dmsref  = Ndj*[30:32].';
thmsref = Ndj*[33:35].';
dssref  = Ndj*[36:38].';
dm1ref  = Ndj*[39:41].';
thm1ref = Ndj*[42:44].';
ThB1ref = Ndj*[45:47].';
dm2ref  = Ndj*[48:50].';
thm2ref = Ndj*[51:53].';
ThB2ref = Ndj*[54:56].';
dm3ref  = Ndj*[57:59].';
thm3ref = Ndj*[60:62].';
ThB3ref = Ndj*[63:65].';
chiref  = Ndj*66;
psiref  = Ndj*67;
epsref  = Ndj*68;
bet1ref = Ndj*69;
bet2ref = Ndj*70;
bet3ref = Ndj*71;

% ================= Foundation-tower joint ===================
% Derivatives are taken with respect to 12 DOFs: thF, dm, thm, and thT, 
% 3 each.
TF0g         = TFromTheta (PB(idofs(1)+[4:6]));
[TFF0,dTFF0] = dTdth      (qB(idofs(1)+[4:6]));
d2TFF0       = d2Tdth2    (qB(idofs(1)+[4:6]),TFF0,dTFF0);
Tm0F         = TFromTheta (PB(idofm(2)+[4:6]));
[Tmm0,dTmm0] = dTdth      (qB(idofm(2)+[4:6]));
d2Tmm0       = d2Tdth2    (qB(idofm(2)+[4:6]),Tmm0,dTmm0);
TTm          = [0 0 1;1 0 0;0 1 0];
TT0g         = TFromTheta (PB(idofs(2)+[4:6]));
[TTT0,dTTT0] = dTdth      (qB(idofs(2)+[4:6]));
d2TTT0       = d2Tdth2    (qB(idofs(2)+[4:6]),TTT0,dTTT0);
TTg_T        = TT0g*TTT0;
TTg_F        = TF0g*TFF0*Tm0F*Tmm0*TTm;  % Necessary for complex step.
TgT          = TTg_T.';
TgF          = Tm0F*Tmm0*TTm*TgT;
Tm0g         = TF0g*TFF0*Tm0F;
Tgm          = TTm*TgT;

ir = 0;  % Displacement constraints.

for jj = 1:3

   jc9 = 9*(jj-1);
   jc3 = 3*(jj-1);

   % ---ThF---
   jref = ThFref(jj);
   dL(ir+[1:3],jref+idofm(2)+[1:3]) = TF0g*dTFF0(:,jc3+[1:3]);

   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(1)+3+kdof) = TF0g*d2TFF0(:,jc9+kd3+[1:3]) ...
                                 * (PB(idofm(2)+[1:3]) + qB(idofm(2)+[1:3]));
   end

   % ---dm---
   jref = dmtref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(1)+3+kdof) = TF0g*dTFF0(:,kd3+jj);
   end

end

ir = ir + 3;  % Rotation constraints.

for jj = 1:3

   jc9 = 9*(jj-1);
   jc3 = 3*(jj-1);

   % ---ThF---
   jref = ThFref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);      

      dS = TF0g*d2TFF0(:,jc9+kd3+[1:3])*TgF;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(1)+3+kdof) = v;

      dS = TF0g*dTFF0(:,jc3+[1:3])*Tm0F*dTmm0(:,kd3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(2)+3+kdof) = v;

      dS = TF0g*dTFF0(:,jc3+[1:3])*Tm0F*Tmm0*TTm*((TT0g*dTTT0(:,kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(2)+3+kdof) = v;

   end

   % ---thm---
   jref = thmtref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);      

      dS = TF0g*dTFF0(:,kd3+[1:3])*Tm0F*dTmm0(:,jc3+[1:3])*TTm*TgT;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(1)+3+kdof) = v;

      dS = Tm0g*d2Tmm0(:,jc9+kd3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(2)+3+kdof) = v;

      dS = Tm0g*dTmm0(:,jc3+[1:3])*TTm*((TT0g*dTTT0(:,kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(2)+3+kdof) = v;

   end

   % ---ThT---
   jref = ThTref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = TF0g*dTFF0(:,kd3+[1:3])*Tm0F*Tmm0*TTm*((TT0g*dTTT0(:,jc3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(1)+3+kdof) = v;

      dS = Tm0g*dTmm0(:,kd3+[1:3])*TTm*((TT0g*dTTT0(:,jc3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(2)+3+kdof) = v;

      dS = TTg_F*((TT0g*d2TTT0(:,jc9+kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(2)+3+kdof) = v;

   end

end

% ================= Tower-nacelle joint ===================
% Derivatives are taken with respect to 13 DOFs: ThT, dm, thm, Thy,
% 3 each, plus chi (yaw angle, scalar).
Tm0T         = TFromTheta (PB(idofm(3)+[4:6]));
[Tmm0,dTmm0] = dTdth      (qB(idofm(3)+[4:6]));
d2Tmm0       = d2Tdth2    (qB(idofm(3)+[4:6]),Tmm0,dTmm0);
Tym          = [0 0 1;1 0 0;0 1 0]*[cy -sy 0;sy cy 0;0 0 1];
dTym         = [0 0 1;1 0 0;0 1 0]*[-sy -cy 0;cy -sy 0;0 0 0];
d2Tym        = [0 0 1;1 0 0;0 1 0]*[-cy sy 0;-sy -cy 0;0 0 0];
Ty0g         = TFromTheta (PB(idofs(3)+[4:6]));
[Tyy0,dTyy0] = dTdth      (qB(idofs(3)+[4:6]));
d2Tyy0       = d2Tdth2    (qB(idofs(3)+[4:6]),Tyy0,dTyy0);
Tyg_Y        = Ty0g*Tyy0;
Tyg_T        = TT0g*TTT0*Tm0T*Tmm0*Tym;  % Necessary for complex step.
Tgy          = Tyg_Y.';
Tm0g         = TT0g*TTT0*Tm0T;
Tgm          = Tym*Tgy;
Tmg          = TT0g*TTT0*Tm0T*Tmm0;      % Necessary for complex step.

ir = ir + 3;  % Displacement constraints.

for jj = 1:3

   jc9 = 9*(jj-1);
   jc3 = 3*(jj-1);

   % ---ThT---
   jref = ThTref(jj);
   dL(ir+[1:3],jref+idofm(3)+[1:3]) = TT0g*dTTT0(:,jc3+[1:3]);
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(2)+3+kdof) = TT0g*d2TTT0(:,jc9+kd3+[1:3]) ...
                               * (PB(idofm(3)+[1:3]) + qB(idofm(3)+[1:3]));
   end

   % ---dm---
   jref = dmyref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(2)+3+kdof) = TT0g*dTTT0(:,kd3+jj);
   end

end

ir = ir + 3;  % Rotation constraints.

for jj = 1:3

   jc9 = 9*(jj-1);
   jc3 = 3*(jj-1);

   % ---ThT---
   jref = ThTref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = TT0g*d2TTT0(:,jc9+kd3+[1:3])*Tm0T*Tmm0*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(2)+3+kdof) = v;

      dS = TT0g*dTTT0(:,jc3+[1:3])*Tm0T*dTmm0(:,kd3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(3)+3+kdof) = v;

      dS = TT0g*dTTT0(:,jc3+[1:3])*Tm0T*Tmm0*Tym*((Ty0g*dTyy0(:,kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(3)+3+kdof) = v;

   end
   dS = TT0g*dTTT0(:,jc3+[1:3])*Tm0T*Tmm0*dTym*Tgy;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(1)) = v;

   % ---thm---
   jref = thmyref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = TT0g*dTTT0(:,kd3+[1:3])*Tm0T*dTmm0(:,jc3+[1:3])*Tym*Tgy;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(2)+3+kdof) = v;

      dS = Tm0g*d2Tmm0(:,jc9+kd3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(3)+3+kdof) = v;

      dS = Tm0g*dTmm0(:,jc3+[1:3])*Tym*((Ty0g*dTyy0(:,kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(3)+3+kdof) = v;

   end
   dS = TT0g*TTT0*Tm0T*dTmm0(:,jc3+[1:3])*dTym*Tgy;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(1)) = v;

   % ---Thy---
   jref = ThYref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = TT0g*dTTT0(:,kd3+[1:3])*Tm0T*Tmm0*Tym*((Ty0g*dTyy0(:,jc3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(2)+3+kdof) = v;

      dS = Tm0g*dTmm0(:,kd3+[1:3])*Tym*((Ty0g*dTyy0(:,jc3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(3)+3+kdof) = v;

      dS = Tyg_T*((Ty0g*d2Tyy0(:,jc9+kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(3)+3+kdof) = v;

   end
   dS = TT0g*TTT0*Tm0T*Tmm0*dTym*((Ty0g*dTyy0(:,jc3+[1:3])).');
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(1)) = v;


end

% ---chi---
jref = chiref;

for kdof = 1:3
   kd3 = 3*(kdof-1);

   dS = TT0g*dTTT0(:,kd3+[1:3])*Tm0T*Tmm0*dTym*Tgy;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofs(2)+3+kdof) = v;

   dS = Tm0g*dTmm0(:,kd3+[1:3])*dTym*Tgy;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofm(3)+3+kdof) = v;

   dS = Tmg*dTym*((Ty0g*dTyy0(:,kd3+[1:3])).');
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofs(3)+3+kdof) = v;

end
dS = Tmg*d2Tym*Tgy;
v = spinToVecU (dS);
dL(ir+[1:3],jref+ijnt(1)) = v;

% ================= Rear bearing ===================
% Derivatives are taken with respect to 13 DOFs: Thy, dm, thm, Thd,
% 3 each, plus psi (azimuth angle, scalar).
Tm0y         = TFromTheta (PB(idofm(4)+[4:6]));
[Tmm0,dTmm0] = dTdth      (qB(idofm(4)+[4:6]));
d2Tmm0       = d2Tdth2    (qB(idofm(4)+[4:6]),Tmm0,dTmm0);
Tdm          = [0 0 -1;-1 0 0;0 1 0]*[ca -sa 0;sa ca 0;0 0 1];
dTdm         = [0 0 -1;-1 0 0;0 1 0]*[-sa -ca 0;ca -sa 0;0 0 0];
d2Tdm        = [0 0 -1;-1 0 0;0 1 0]*[-ca sa 0;-sa -ca 0;0 0 0];
Td0g         = TFromTheta (PB(idofs(4)+[4:6]));
[Tdd0,dTdd0] = dTdth      (qB(idofs(4)+[4:6]));
d2Tdd0       = d2Tdth2    (qB(idofs(4)+[4:6]),Tdd0,dTdd0);
Tdg_D        = Td0g*Tdd0;
Tdg_Y        = Ty0g*Tyy0*Tm0y*Tmm0*Tdm;  % Necessary for complex step.
Tgd          = Tdg_D.';
Tm0g         = Ty0g*Tyy0*Tm0y;
Tgm          = Tdm*Tgd;
Tmg          = Ty0g*Tyy0*Tm0y*Tmm0;      % Necessary for complex step.

ir = ir + 3;  % Displacement constraints.

for jj = 1:3

   jc9 = 9*(jj-1);
   jc3 = 3*(jj-1);

   % ---Thy---
   jref = ThYref(jj);
   dL(ir+[1:3],jref+idofm(4)+[1:3]) = Ty0g*dTyy0(:,jc3+[1:3]);
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(3)+3+kdof) = Ty0g*d2Tyy0(:,jc9+kd3+[1:3]) ...
                                 * (PB(idofm(4)+[1:3]) + qB(idofm(4)+[1:3]));
   end

   % ---dm---
   jref = dmdref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(3)+3+kdof) = Ty0g*dTyy0(:,kd3+jj);
   end

end

ir = ir + 3;  % Rotation constraints.

for jj = 1:3

   jc9 = 9*(jj-1);
   jc3 = 3*(jj-1);

   % ---Thy---
   jref = ThYref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = Ty0g*d2Tyy0(:,jc9+kd3+[1:3])*Tm0y*Tmm0*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(3)+3+kdof) = v;

      dS = Ty0g*dTyy0(:,jc3+[1:3])*Tm0y*dTmm0(:,kd3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(4)+3+kdof) = v;

      dS = Ty0g*dTyy0(:,jc3+[1:3])*Tm0y*Tmm0*Tdm*((Td0g*dTdd0(:,kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

   end
   dS = Ty0g*dTyy0(:,jc3+[1:3])*Tm0y*Tmm0*dTdm*Tgd;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(2)) = v;

   % ---thm---
   jref = thmdref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = Ty0g*dTyy0(:,kd3+[1:3])*Tm0y*dTmm0(:,jc3+[1:3])*Tdm*Tgd;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(3)+3+kdof) = v;

      dS = Ty0g*Tyy0*Tm0y*d2Tmm0(:,jc9+kd3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(4)+3+kdof) = v;

      dS = Ty0g*Tyy0*Tm0y*dTmm0(:,jc3+[1:3])*Tdm*((Td0g*dTdd0(:,kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

   end
   dS = Ty0g*Tyy0*Tm0y*dTmm0(:,jc3+[1:3])*dTdm*Tgd;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(2)) = v;

   % ---Thd---
   jref = ThDref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = Ty0g*dTyy0(:,kd3+[1:3])*Tm0y*Tmm0*Tdm*((Td0g*dTdd0(:,jc3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(3)+3+kdof) = v;

      dS = Tm0g*dTmm0(:,kd3+[1:3])*Tdm*((Td0g*dTdd0(:,jc3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(4)+3+kdof) = v;

      dS = Tdg_Y*((Td0g*d2Tdd0(:,jc9+kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

   end
   dS = Tmg*dTdm*((Td0g*dTdd0(:,jc3+[1:3])).');
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(2)) = v;

end

% ---psi---
jref = psiref;
for kdof = 1:3
   kd3 = 3*(kdof-1);

   dS = Ty0g*dTyy0(:,kd3+[1:3])*Tm0y*Tmm0*dTdm*Tgd;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofs(3)+3+kdof) = v;

   dS = Tm0g*dTmm0(:,kd3+[1:3])*dTdm*Tgd;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofm(4)+3+kdof) = v;

   dS = Tmg*dTdm*((Td0g*dTdd0(:,kd3+[1:3])).');
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

end
dS = Tmg*d2Tdm*Tgd;
v = spinToVecU (dS);
dL(ir+[1:3],jref+ijnt(2)) = v;

% ================= Front bearing ===================
% Derivatives are taken with respect to 16 DOFs: Thy, dm, thm, Thd,
% ds, 3 each, and eps, scalar.
Tm0y         = TFromTheta (PB(idofm(5)+[4:6]));
[Tmm0,dTmm0] = dTdth      (qB(idofm(5)+[4:6]));
d2Tmm0       = d2Tdth2    (qB(idofm(5)+[4:6]),Tmm0,dTmm0);

ir = ir + 3;  % Only displacement constraints, in this case.

for jj = 1:3

   jc9 = 9*(jj-1);
   jc3 = 3*(jj-1);

   % ---Thy---
   jref = ThYref(jj);
   dL(ir+[1:3],jref+idofm(5)+[1:3]) = Ty0g*dTyy0(:,jc3+[1:3]);
   dL(ir+[1:3],jref+ijnt(3)) = Ty0g*dTyy0(:,jc3+[1:3])*Tm0y*Tmm0*[1;0;0];
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(3)+3+kdof) = Ty0g*d2Tyy0(:,jc9+kd3+[1:3]) ...
          *(PB(idofm(5)+[1:3]) + qB(idofm(5)+[1:3]) + Tm0y*Tmm0*[fbx;0;0]);
      dL(ir+[1:3],jref+idofm(5)+3+kdof) = Ty0g*dTyy0(:,jc3+[1:3]) ...
                                        * Tm0y*dTmm0(:,kd3+[1:3])*[fbx;0;0];
   end

   % ---dm---
   jref = dmsref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(3)+3+kdof) = Ty0g*dTyy0(:,kd3+jj);
   end

   % ---thm---
   jref = thmsref(jj);
   dL(ir+[1:3],jref+ijnt(3)) = Tyg_Y*Tm0y*dTmm0(:,jc3+[1:3])*[1;0;0];
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(3)+3+kdof) = Ty0g*dTyy0(:,kd3+[1:3]) ...
                                        * Tm0y*dTmm0(:,jc3+[1:3])*[fbx;0;0];
      dL(ir+[1:3],jref+idofm(5)+3+kdof) = Tyg_Y*Tm0y*d2Tmm0(:,jc9+kd3+[1:3])*[fbx;0;0];
   end

   % ---Thd---
   jref = ThDref(jj);
   dL(ir+[1:3],jref+idofs(5)+[1:3]) = -Td0g*dTdd0(:,jc3+[1:3]);
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = -Td0g*d2Tdd0(:,jc9+kd3+[1:3]) ...
                               * (PB(idofs(5)+[1:3]) + qB(idofs(5)+[1:3]));
   end

   % ---ds---
   jref = dssref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = -Td0g*dTdd0(:,kd3+jj);
   end

end

% ---eps---
jref = epsref;
for kdof = 1:3
   kd3 = 3*(kdof-1);
   dL(ir+[1:3],jref+idofs(3)+3+kdof) = Ty0g*dTyy0(:,kd3+[1:3])*Tm0y*Tmm0*[1;0;0];
   dL(ir+[1:3],jref+idofm(5)+3+kdof) = Tyg_Y*Tm0y*dTmm0(:,kd3+[1:3])*[1;0;0];
end

% ================= Blade 1 ===================
% Derivatives are taken with respect to 13 DOFs: Thd, dm, thm, Thp,
% 3 each, and beta, scalar.
Tm0d         = TFromTheta (PB(idofm(6)+[4:6]));
[Tmm0,dTmm0] = dTdth      (qB(idofm(6)+[4:6]));
d2Tmm0       = d2Tdth2    (qB(idofm(6)+[4:6]),Tmm0,dTmm0);
Tpm          = Tbh*[1 0 0;0 cb1 sb1;0 -sb1 cb1];
dTpm         = Tbh*[0 0 0;0 -sb1 cb1;0 -cb1 -sb1];
d2Tpm        = Tbh*[0 0 0;0 -cb1 -sb1;0 sb1 -cb1];
Tp0g         = TFromTheta (PB(idofs(6)+[4:6]));
[Tpp0,dTpp0] = dTdth      (qB(idofs(6)+[4:6]));
d2Tpp0       = d2Tdth2    (qB(idofs(6)+[4:6]),Tpp0,dTpp0);
Tpg_P        = Tp0g*Tpp0;
Tpg_D        = Td0g*Tdd0*Tm0d*Tmm0*Tpm;  % Necessary for complex step.
Tgp          = Tpg_P.';
Tm0g         = Td0g*Tdd0*Tm0d;
Tgm          = Tpm*Tgp;
Tmg          = Td0g*Tdd0*Tm0d*Tmm0;      % Necessary for complex step.

ir = ir + 3;  % Displacement constraints.

for jj = 1:3

   jc9 = 9*(jj-1);
   jc3 = 3*(jj-1);

   % ---Thd---
   jref = ThDref(jj);
   dL(ir+[1:3],jref+idofm(6)+[1:3]) = Td0g*dTdd0(:,jc3+[1:3]);
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = Td0g*d2Tdd0(:,jc9+kd3+[1:3]) ...
                               * (PB(idofm(6)+[1:3]) + qB(idofm(6)+[1:3]));
   end

   % ---dm---
   jref = dm1ref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = Td0g*dTdd0(:,kd3+jj);
   end

end

ir = ir + 3;  % Rotation constraints.

for jj = 1:3

   jc9 = 9*(jj-1);
   jc3 = 3*(jj-1);

   % ---Thd---
   jref = ThDref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = Td0g*d2Tdd0(:,jc9+kd3+[1:3])*Tm0d*Tmm0*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

      dS = Td0g*dTdd0(:,jc3+[1:3])*Tm0d*dTmm0(:,kd3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(6)+3+kdof) = v;

      dS = Td0g*dTdd0(:,jc3+[1:3])*Tm0d*Tmm0*Tpm*((Tp0g*dTpp0(:,kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(6)+3+kdof) = v;

   end
   dS = Td0g*dTdd0(:,jc3+[1:3])*Tm0d*Tmm0*dTpm*Tgp;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(4)) = v;

   % ---thm---
   jref = thm1ref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = Td0g*dTdd0(:,kd3+[1:3])*Tm0d*dTmm0(:,jc3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

      dS = Tm0g*d2Tmm0(:,jc9+kd3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(6)+3+kdof) = v;

      dS = Tm0g*dTmm0(:,jc3+[1:3])*Tpm*((Tp0g*dTpp0(:,kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(6)+3+kdof) = v;

   end
   dS = Tm0g*dTmm0(:,jc3+[1:3])*dTpm*Tgp;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(4)) = v;

   % ---ThB1---
   jref = ThB1ref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = Td0g*dTdd0(:,kd3+[1:3])*Tm0d*Tmm0*Tpm*((Tp0g*dTpp0(:,jc3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

      dS = Tm0g*dTmm0(:,kd3+[1:3])*Tpm*((Tp0g*dTpp0(:,jc3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(6)+3+kdof) = v;

      dS = Tmg*Tpm*((Tp0g*d2Tpp0(:,jc9+kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(6)+3+kdof) = v;

   end
   dS = Tmg*dTpm*((Tp0g*dTpp0(:,jc3+[1:3])).');
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(4)) = v;

end

% ---beta1---
jref = bet1ref;
for kdof = 1:3
   kd3 = 3*(kdof-1);

   dS = Td0g*dTdd0(:,kd3+[1:3])*Tm0d*Tmm0*dTpm*Tgp;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

   dS = Tm0g*dTmm0(:,kd3+[1:3])*dTpm*Tgp;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofm(6)+3+kdof) = v;

   dS = Tmg*dTpm*((Tp0g*dTpp0(:,kd3+[1:3])).');
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofs(6)+3+kdof) = v;

end
dS = Tmg*d2Tpm*Tgp;
v = spinToVecU (dS);
dL(ir+[1:3],jref+ijnt(4)) = v;

% ================= Blade 2 ===================
% Derivatives are taken with respect to 13 DOFs: Thd, dm, thm, Thp,
% 3 each, and beta, scalar.
Tm0d         = TFromTheta (PB(idofm(7)+[4:6]));
[Tmm0,dTmm0] = dTdth      (qB(idofm(7)+[4:6]));
d2Tmm0       = d2Tdth2    (qB(idofm(7)+[4:6]),Tmm0,dTmm0);
Tpm          = Tbh*[1 0 0;0 cb2 sb2;0 -sb2 cb2];
dTpm         = Tbh*[0 0 0;0 -sb2 cb2;0 -cb2 -sb2];
d2Tpm        = Tbh*[0 0 0;0 -cb2 -sb2;0 sb2 -cb2];
Tp0g         = TFromTheta (PB(idofs(7)+[4:6]));
[Tpp0,dTpp0] = dTdth      (qB(idofs(7)+[4:6]));
d2Tpp0       = d2Tdth2    (qB(idofs(7)+[4:6]),Tpp0,dTpp0);
Tpg_P        = Tp0g*Tpp0;
Tpg_D        = Td0g*Tdd0*Tm0d*Tmm0*Tpm;  % Necessary for complex step.
Tgp          = Tpg_P.';
Tm0g         = Td0g*Tdd0*Tm0d;
Tgm          = Tpm*Tgp;
Tmg          = Td0g*Tdd0*Tm0d*Tmm0;      % Necessary for complex step.

ir = ir + 3;  % Displacement constraints.

for jj = 1:3

   jc9 = 9*(jj-1);
   jc3 = 3*(jj-1);

   % ---Thd---
   jref = ThDref(jj);
   dL(ir+[1:3],jref+idofm(7)+[1:3]) = Td0g*dTdd0(:,jc3+[1:3]);
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = Td0g*d2Tdd0(:,jc9+kd3+[1:3]) ...
                               * (PB(idofm(7)+[1:3]) + qB(idofm(7)+[1:3]));
   end

   % ---dm---
   jref = dm2ref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = Td0g*dTdd0(:,kd3+jj);
   end

end

ir = ir + 3;  % Rotation constraints.

for jj = 1:3

   jc9 = 9*(jj-1);
   jc3 = 3*(jj-1);

   % ---Thd---
   jref = ThDref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = Td0g*d2Tdd0(:,jc9+kd3+[1:3])*Tm0d*Tmm0*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

      dS = Td0g*dTdd0(:,jc3+[1:3])*Tm0d*dTmm0(:,kd3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(7)+3+kdof) = v;

      dS = Td0g*dTdd0(:,jc3+[1:3])*Tm0d*Tmm0*Tpm*((Tp0g*dTpp0(:,kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(7)+3+kdof) = v;

   end
   dS = Td0g*dTdd0(:,jc3+[1:3])*Tm0d*Tmm0*dTpm*Tgp;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(5)) = v;

   % ---thm---
   jref = thm2ref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = Td0g*dTdd0(:,kd3+[1:3])*Tm0d*dTmm0(:,jc3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

      dS = Tm0g*d2Tmm0(:,jc9+kd3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(7)+3+kdof) = v;

      dS = Tm0g*dTmm0(:,jc3+[1:3])*Tpm*((Tp0g*dTpp0(:,kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(7)+3+kdof) = v;

   end
   dS = Tm0g*dTmm0(:,jc3+[1:3])*dTpm*Tgp;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(5)) = v;

   % ---ThB2---
   jref = ThB2ref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = Td0g*dTdd0(:,kd3+[1:3])*Tm0d*Tmm0*Tpm*((Tp0g*dTpp0(:,jc3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

      dS = Tm0g*dTmm0(:,kd3+[1:3])*Tpm*((Tp0g*dTpp0(:,jc3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(7)+3+kdof) = v;

      dS = Tmg*Tpm*((Tp0g*d2Tpp0(:,jc9+kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(7)+3+kdof) = v;

   end
   dS = Tmg*dTpm*((Tp0g*dTpp0(:,jc3+[1:3])).');
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(5)) = v;

end

% ---beta2---
jref = bet2ref;
for kdof = 1:3
   kd3 = 3*(kdof-1);

   dS = Td0g*dTdd0(:,kd3+[1:3])*Tm0d*Tmm0*dTpm*Tgp;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

   dS = Tm0g*dTmm0(:,kd3+[1:3])*dTpm*Tgp;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofm(7)+3+kdof) = v;

   dS = Tmg*dTpm*((Tp0g*dTpp0(:,kd3+[1:3])).');
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofs(7)+3+kdof) = v;

end
dS = Tmg*d2Tpm*Tgp;
v = spinToVecU (dS);
dL(ir+[1:3],jref+ijnt(5)) = v;

% ================= Blade 3 ===================
% Derivatives are taken with respect to 13 DOFs: Thd, dm, thm, Thp,
% 3 each, and beta, scalar.
Tm0d         = TFromTheta (PB(idofm(8)+[4:6]));
[Tmm0,dTmm0] = dTdth      (qB(idofm(8)+[4:6]));
d2Tmm0       = d2Tdth2    (qB(idofm(8)+[4:6]),Tmm0,dTmm0);
Tpm          = Tbh*[1 0 0;0 cb3 sb3;0 -sb3 cb3];
dTpm         = Tbh*[0 0 0;0 -sb3 cb3;0 -cb3 -sb3];
d2Tpm        = Tbh*[0 0 0;0 -cb3 -sb3;0 sb3 -cb3];
Tp0g         = TFromTheta (PB(idofs(8)+[4:6]));
[Tpp0,dTpp0] = dTdth      (qB(idofs(8)+[4:6]));
d2Tpp0       = d2Tdth2    (qB(idofs(8)+[4:6]),Tpp0,dTpp0);
Tpg_P        = Tp0g*Tpp0;
Tpg_D        = Td0g*Tdd0*Tm0d*Tmm0*Tpm;  % Necessary for complex step.
Tgp          = Tpg_P.';
Tm0g         = Td0g*Tdd0*Tm0d;
Tgm          = Tpm*Tgp;
Tmg          = Td0g*Tdd0*Tm0d*Tmm0;      % Necessary for complex step.

ir = ir + 3;  % Displacement constraints.

for jj = 1:3

   jc9 = 9*(jj-1);
   jc3 = 3*(jj-1);

   % ---Thd---
   jref = ThDref(jj);
   dL(ir+[1:3],jref+idofm(8)+[1:3]) = Td0g*dTdd0(:,jc3+[1:3]);
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = Td0g*d2Tdd0(:,jc9+kd3+[1:3]) ...
                               * (PB(idofm(8)+[1:3]) + qB(idofm(8)+[1:3]));
   end

   % ---dm---
   jref = dm3ref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = Td0g*dTdd0(:,kd3+jj);
   end

end

ir = ir + 3;  % Rotation constraints.

for jj = 1:3

   jc9 = 9*(jj-1);
   jc3 = 3*(jj-1);

   % ---Thd---
   jref = ThDref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = Td0g*d2Tdd0(:,jc9+kd3+[1:3])*Tm0d*Tmm0*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

      dS = Td0g*dTdd0(:,jc3+[1:3])*Tm0d*dTmm0(:,kd3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(8)+3+kdof) = v;

      dS = Td0g*dTdd0(:,jc3+[1:3])*Tm0d*Tmm0*Tpm*((Tp0g*dTpp0(:,kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(8)+3+kdof) = v;

   end
   dS = Td0g*dTdd0(:,jc3+[1:3])*Tm0d*Tmm0*dTpm*Tgp;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(6)) = v;

   % ---thm---
   jref = thm3ref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = Td0g*dTdd0(:,kd3+[1:3])*Tm0d*dTmm0(:,jc3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

      dS = Tm0g*d2Tmm0(:,jc9+kd3+[1:3])*Tgm;
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(8)+3+kdof) = v;

      dS = Tm0g*dTmm0(:,jc3+[1:3])*Tpm*((Tp0g*dTpp0(:,kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(8)+3+kdof) = v;

   end
   dS = Tm0g*dTmm0(:,jc3+[1:3])*dTpm*Tgp;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(6)) = v;

   % ---ThB3---
   jref = ThB3ref(jj);
   for kdof = 1:3
      kd3 = 3*(kdof-1);

      dS = Td0g*dTdd0(:,kd3+[1:3])*Tm0d*Tmm0*Tpm*((Tp0g*dTpp0(:,jc3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

      dS = Tm0g*dTmm0(:,kd3+[1:3])*Tpm*((Tp0g*dTpp0(:,jc3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofm(8)+3+kdof) = v;

      dS = Tmg*Tpm*((Tp0g*d2Tpp0(:,jc9+kd3+[1:3])).');
      v = spinToVecU (dS);
      dL(ir+[1:3],jref+idofs(8)+3+kdof) = v;

   end
   dS = Tmg*dTpm*((Tp0g*dTpp0(:,jc3+[1:3])).');
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+ijnt(6)) = v;

end

% ---beta3---
jref = bet3ref;
for kdof = 1:3
   kd3 = 3*(kdof-1);

   dS = Td0g*dTdd0(:,kd3+[1:3])*Tm0d*Tmm0*dTpm*Tgp;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofs(4)+3+kdof) = v;

   dS = Tm0g*dTmm0(:,kd3+[1:3])*dTpm*Tgp;
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofm(8)+3+kdof) = v;

   dS = Tmg*dTpm*((Tp0g*dTpp0(:,kd3+[1:3])).');
   v = spinToVecU (dS);
   dL(ir+[1:3],jref+idofs(8)+3+kdof) = v;

end
dS = Tmg*d2Tpm*Tgp;
v = spinToVecU (dS);
dL(ir+[1:3],jref+ijnt(6)) = v;

