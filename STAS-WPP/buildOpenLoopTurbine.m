function [L,R,yout,A,B,C,D] =                                              ...
               buildOpenLoopTurbine (linFlag,x,u,s,a,epar,ppar,ypar,       ...
                                     grav,P,shape0,mdamp0,                 ...
                                     Tas,Try,ch,Lel,foilwt,aoaz,aoast,     ...
                                     xas,yas,Psi0,igen,ipit,iyaw)
%
% Build the coupled aeroelastic, electrical, and actuator models, such that
% inputs are the windspeed and actuator control commands. 
%
%   States:              y vector:             u vector:
% ----------------------- Structure --------------------------
%   eta        N         q         Ndj            F       Ndj   (Env)
%   deta/dt    N         dq/dt     Ndj         d2eta/dt2  Neta  (ext. solution)
%                        d2q/dt2   Ndj
%                        xng     3*Nnod
%                        vng     6*Nnod
%                        F         Ndj
%
% ---------------------- Aerodynamic -------------------------
%   ad         1         (q)       1:24           Vg        3   (Env)
%   a1,a2     2:3        (dq/dt)  25:48
%   Vih z,t   4:5        xng1,2   49:54
%   Vi z,t    6:7        vng1,2   55:60
%                        wg       61:63
%                        wa       64:66
%                        aq        67
%                        Cl,Cd,Cm 68:70
%                        Fl,Fd,M  71:73
%                        Fa       74:79
%                        Fp       80:85
%                        Fr       86:91
%                        Fzts     92:94
%                        Vg       95:97
%                        Ua       98:100
%                        Umag      101
%                        Ur      102:104
%                        Uzts    105:107
%                        Vzts    108:110
%                        Viq     111:112
%                        Viy     113:114
%                        Vixyz   115:117
%                        Wmag      118
%                        xeg     119:121
%                        xhg     122:124
%                        xnr1,2  125:130  
%                        xer     131:133
%                        r         134
%                        Lp        135
%                        z         136
%                        f         137  
% (Repeat the above x,y,u consecutively for each blade element.)
%                        Azi       (1)  (effective rotor azimuth)
%                        Waero     (1)  (aero rotor speed)
%                        Dp        (3)
%
% ------------------------ Actuators --------------------------
% (Repeat for each blade and yaw system.)
%   ba         1         b,dbdt    1,2  (Turb)     bhat      1   (Cont)
%   dbadt      2         Ta         3
%
% ----------------------- Electrical --------------------------
%   igd,q     1,2        wg         1   (Turb)     we        1   (Grid)
%   Vdc        3         Tg         2              th_e      2   (Grid)
%   ipd,q     4,5        isd,q     3,4             vsd,q    3,4  (Grid)
%   imgd,q    6,7                                  ihgd,q   5,6  (Cont)
%   Psig      8,9                                  Vhdc      7   (Cont)
%   wemg      10                                   Qh        8   (Cont)
%   th_m      11
%   vmsd,q   12,13
%   Psie      14
%   impd,q   15,16
%   imsd,q   17,18
%   vmsd,q   19,20
%   Vmdc      21
%   PsiDC     22
%   PsiQ      23
%   PsiP     24,25
%
% Version:        Changes:
% --------        -------------
% 12.01.2019      Original code.
% 22.01.2020      Minor updates to accommodate revised aeroelastic.m,
%                 including updated structural equations-of-motion.
%
% Version:        Verification:
% --------        -------------
% 12.01.2019      Most terms in the A matrix have been verified by complex
%                 step.  There is a known issue with the nacelle yaw position,
%                 here the linear and nonlinear equations don't agree.
% 22.01.2020      "A" matrix derivatives verified using complex step on L,R.
%                 The above issue with nacelle yaw was fixed.
%
% Inputs:
% -------
% Linflag         : Set to 1 for linearized equations.
% x               :  N    eta   (-)     Modal displacements.
%                    N    deta  (1/s)   Rate of change of modal displacements.   
%                   ---
%                    1    ad    (rad)   Dynamic angle-of-attack.
%    Each blade     2,3   a1,a2         Intermediate aoa variables.
%    element        4,5   Vih           Intermediate ind. velocity variables.
%                   6,7   Vi    (m/s)   Induced velocity.
%                   ---
%    Each actuator   1    ba    (rad)   Actuator angle.
%                    2    dbadt (rad/s) Actuator speed.
%                   ---
%                   1,2   igd,q (A)     Generator d,q axis currents.
%                    3    Vdc   (V)     DC link voltage.
%                   4,5   ipd,q (A)     Currents at the network-side converter.
%                   6,7   img   (A)     Measured generator currents.
%                   8,9   Psig  (As)    Integrated generator currents.
%                   10    wemg  (rad/s) Measured generator speed.
%                   11    th_m  (rad)   Measured electrical angle.
%                  12,13  vms   (V)     Measured transformer terminal voltage.
%                   14    Psie  (rad s) Integrated electrical angle.
%                  15,16  imp   (A)     Measured primary terminal current.
%                  17,18  ims   (A)     Measured secondary terminal current.
%                  19,20  vms   (V)     Measured secondary terminal voltage.
%                   21    Vmdc  (V)     Measured DC link voltage.
%                   22    PsiDC (Vs)    Integrated DC link voltage error.
%                   23    PsiQ  (VAs)   Integrated reactive power error.
%                  24,25  PsiP  (Ws)    Integrated current error.
%
% u               : Ndj   F     (N,Nm)  Applied nodal forces and moments.
%                   Neta  etadd (1/s^2) Rate-of-rate-of-change of...
%                   ---
%    Each bld. el.   3    Vg    (m/s)   Global V wind for each blade element.
%                   ---
%    Each actuator   1    bhat  (rad)   Actuator command for each joint.
%                   ---
%                    1    we    (rad/s) Electrical speed of network.
%                    2    th_e  (rad)   Electrical angle of network.
%                   3,4   vsd,q (V)     Voltage at network-side transformer terminals.
%                   5,6   ihg   (A)     Generator current command.
%                    7    Vhdc  (V)     Commanded DC link voltage.
%                    8    Qh    (VA)    Reactive power command.
%                    
% epar            : Parameters for the electrical system, see buildTurbineElectric.m
% ppar            : Parameters for the pitch actuators, see buildDrive.m
% ypar            : Parameters for the yaw actuator, see buildDrive.m
%
% Outputs:
% --------
% dxdt            : Rate of change of states.
% yout            : The full y vector.
% A,B,C,D         : Linearized state matrices.

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

Nnod = Ndof/6;
Ndj = size(P,1);
Neta = size(shape0,2);
Neb = a.Neb;
Nae = a.Neb*a.Nb;

Nx  = size(x,1);
Nxs = 2*Neta;
Nxa = size(Psi0,2);
Nx1 = Nxs + Nxa;
Nxp = 6;
Nxy = 2;
Nxe = 25;

Nu  = size(u,1);
Nus = Ndj + Neta;
Nua = 3*Nae;
Nup = a.Nb;
Nuy = 1;
Nue = 8;

ixs = 0;
ixa = ixs + Nxs;
ixp = ixa + Nxa;
ixy = ixp + Nxp;
ixe = ixy + Nxy;

ius = 0;
iua = ius + Nus;
iup = iua + Nua;
iuy = iup + Nup;
iue = iuy + Nuy;

% Extract q's from x's.
slv = slaveDOFs (idofs);
vec = [1:Ndj].';
[jnk,ret,jnk2] = partitionMatrix (vec,slv,[]);
xs   = x(1:Nxs);
eta  = xs(1:Neta);
etad = xs(Neta+[1:Neta]);
zro  = zeros(Neta,1);     % Not needed for nonlinear.

[qh,qhd,jnk1,jnk2] = buildTetaqh (eta,etad,zro,shape0);

[q,dqdt,jnk2,jnk3,jnk4,jnk5] = buildTqhq (0,s,qh,qhd,jnk1,P,ret,slv);

% Electrical model.
qguess = zeros(size(slv,1),1);
[q,dqdt] = qFromx (x,s,P,shape0,ret,slv,qguess);
qgen = q(igen);
dqgen = dqdt(igen);
Pgen = P(igen);
[GenSp,ddyw] = genSpeed (qgen,dqgen,Pgen);
xe = x(ixe+[1:Nxe]);
yein = [0.5*epar(1)*GenSp;u(iue+[1:8])];
[dxedt,yeout,aae,bbye,cce,ddye] = buildTurbineElectric (linFlag,xe,yein,epar);

% Pitch actuators.
dxpdt = zeros(Nxp,1);
ypout = zeros(a.Nb,1);
for ib = 1:a.Nb
   ir2 = 2*(ib-1);
   ir4 = 4*(ib-1);
   indx = ir2 + [1:2].';
   indy = ir4 + [1:4].';
   xp = x(ixp+ir2+[1:2]);
   ypin = [u(iup+ib);q(ipit(ib));dqdt(ipit(ib))];
   [dxpdt(indx),ypout(ib),aap(indx,indx),bbyp(indx,indy), ...
    ccp(indy,indx),ddyp(indy,indy)] = buildDrive (linFlag,xp,ypin,ppar);
end

% Yaw actuator.
xy = x(ixy+[1:2]);
yyin = [u(iuy+1);q(iyaw);dqdt(iyaw)];
[dxydt,yyout,aay,bbyy,ccy,ddyy] = buildDrive (linFlag,xy,yyin,ypar);

% Aeroelastic model.
x1  = x(1:Nx1);

% Apply the generator torque.
Fext = zeros(Ndj,1);
Tmm0 = TFromTheta (qgen(10:12));
Tm0y = TFromTheta (Pgen(10:12));
Tnn0 = TFromTheta (qgen(22:24));
Tn0d = TFromTheta (Pgen(22:24));
Fext(igen(10:12)) = Tm0y*Tmm0*[-yeout(1);0;0];
Fext(igen(22:24)) = Tn0d*Tnn0*[yeout(1);0;0];

% Apply the pitch actuator torques.
for ib = 1:3
   Tmm0 = TFromTheta (q(idofm(5+ib)+[4:6]));
   Tm0d = TFromTheta (P(idofm(5+ib)+[4:6]));
   Fext(idofm(5+ib)+[4:6]) = Tm0d*Tmm0*[ypout(ib);0;0];
   Fext(idofs(5+ib)+4) = -ypout(ib);
end

% Apply the yaw actuator torque.
Tmm0 = TFromTheta (q(idofm(3)+[4:6]));
Tm0t = TFromTheta (P(idofm(3)+[4:6]));
Fext(idofm(3)+[4:6]) = Tm0t*Tmm0*[-yyout;0;0];
Fext(idofs(3)+6) = yyout;

u1        = u(1:Nus+Nua);
u1(1:Ndj) = u1(1:Ndj) + Fext;

[Lmat,Rvec,yae,jA,jBu,jBy,jC,jDu,jDy] =   ...
   aeroelastic (0,s,a,x1,u1,P,            ...
                shape0,mdamp0,grav,       ...
                Tas,Try,ch,Lel,foilwt,    ...
                aoaz,aoast,xas,yas,Psi0);
Nyal = size(jC,1);

% Form the nonlinear outputs.
indy = [[1:3*Ndj] 3*Ndj+6*Nae+[1:Ndj]].';

Nxpye = Nxp + Nxy + Nxe;
L = [Lmat sparse(Nx1,Nxpye);sparse(Nxpye,Nx1) speye(Nxpye)];
R = [Rvec;dxpdt;dxydt;dxedt];

yout = [yae(indy);ypout;yyout;yeout];
Nyout = size(yout,1);

if (linFlag == 1)

   % Generate the linearized aeroelastic matrices.
   q = yae(1:Ndj);
   dq = yae(Ndj+[1:Ndj]);
   d2q = yae(2*Ndj+[1:Ndj]);
   F = yae(3*Ndj+6*Nae+[1:Ndj]);

   [Lmat,Rvec,yae,aa1,bbu1,bby1,cc1,ddu1,ddy1] = ...
         aeroelastic (1,s,a,x1,u1,P,             ...
                      shape0,mdamp0,grav,        ...
                      Tas,Try,ch,Lel,foilwt,     ...
                      aoaz,aoast,xas,yas,Psi0);

   Ny1 = size(cc1,1);
   Nyp = size(ccp,1) - Nup;  % Going to partition out u variables.
   Nyy = size(ccy,1) - Nuy;
   Nye = size(yeout,1) + size(yein,1) - Nue;
   Ny = Ny1 + Nyp + Nyy + Nye;

   iy1 = 0;
   iyp = iy1 + Ny1;
   iyy = iyp + Nyp;
   iye = iyy + Nyy;

   % Fill the un-linked matrices. 
   aa  = spalloc (Nx,Nx,round(0.3*Nx*Nx));
   bbu = spalloc (Nx,Nu,round(0.15*Nx*Nu));
   bby = spalloc (Nx,Ny,round(0.02*Nx*Ny));
   cc  = spalloc (Ny,Nx,round(0.08*Ny*Nx));
   ddu = spalloc (Ny,Nu,round(0.08*Ny*Nu));
   ddy = spalloc (Ny,Ny,round(1e-5*Ny*Ny));

   indx = ixs  + [1:Nx1].';
   indu = ius  + [1:Nus+Nua].';
   indy = iy1  + [1:Ny1].';
   aa(indx,indx)  = aa1;
   bbu(indx,indu) = bbu1;
   bby(indx,indy) = bby1;
   cc(indy,indx)  = cc1;
   ddu(indy,indu) = ddu1;
   ddy(indy,indy) = ddy1;

   indx = ixp  + [1:Nxp].';
   indu = iup  + [1:Nup].';
   indy = iyp  + [1:Nyp].';
   ju = [1 5 9].';
   jy = [2 3 4 6 7 8 10 11 12].';
   aa(indx,indx)  = aap;
   bbu(indx,indu) = bbyp(:,ju);
   bby(indx,indy) = bbyp(:,jy);
   cc(indy,indx)  = ccp(jy,:);
   ddu(indy,indu) = ddyp(jy,ju);
   ddy(indy,indy) = ddyp(jy,jy);

   indx = ixy  + [1:Nxy].';
   indu = iuy  + [1:Nuy].';
   indy = iyy  + [1:Nyy].';
   ju = 1;
   jy = [2 3 4].';
   aa(indx,indx)  = aay;
   bbu(indx,indu) = bbyy(:,ju);
   bby(indx,indy) = bbyy(:,jy);
   cc(indy,indx)  = ccy(jy,:);
   ddu(indy,indu) = ddyy(jy,ju);
   ddy(indy,indy) = ddyy(jy,jy);

   % The indexing in buildTurbineElectric is handled a bit differently,
   % distinguishing between inputs and outputs.
   indx = ixe  + [1:Nxe].';
   indu = iue  + [1:Nue].';   % The turbine indices for the u inputs.
   indyi = iye + 1;           % The turbine index of electrical input wg.
   indyo = iye + [2:4].';     % The turbine indices of electrical outputs Tg, isd,q.
   ju = [2:9].';              % The elec u (input) indices for the turbine u inputs.
   jyi = 1;                   % The electric u (input) index for wg.
   jyo = [1:3].';             % The electric y (output) indices for Tg, isd,q. 
   aa(indx,indx)    = aae;
   bbu(indx,indu)   = bbye(:,ju);
   bby(indx,indyi)  = bbye(:,jyi);
   cc(indyo,indx)   = cce(jyo,:);
   ddu(indyo,indu)  = ddye(jyo,ju);
   ddy(indyo,indyi) = ddye(jyo,jyi);

   % Link the turbine blade pitch and pitch rate to the actuators.
   ir = iyp + [1 4 7 2 5 8].';
   ic = iy1 + [ipit;Ndj+ipit];
   ddy(ir,ic) = speye(6);

   % Link the pitch actuator torque to the blade root nodes.
   for ib = 1:3
      mdofs = idofm(5+ib)+[4:6];
      [Tmm0,dTmm0] = dTdth (q(mdofs));
      Tm0d = TFromTheta (P(mdofs));
      ir = 3*Ndj + 9*Nnod + mdofs;
      ic = iyp + 3*ib;
      ddy(ir,ic) = Tm0d*Tmm0(:,1);
      ic = mdofs;
      for ith = 1:3
         ic3 = 3*(ith-1);
         ddy(ir,ic(ith)) = ddy(ir,ic(ith)) ...
                         + Tm0d*dTmm0(:,ic3+1)*ypout(ib);
      end
      ir = 3*Ndj + 9*Nnod + idofs(5+ib)+4;
      ic = iyp + 3*ib;
      ddy(ir,ic) = -1;
   end

   % Link the turbine yaw angle and yaw rate to the actuator.
   ir = iyy + [1;2];
   ic = iy1 + [iyaw;Ndj+iyaw];
   ddy(ir,ic) = speye(2);

   % Link the yaw actuator torque to the tower and nacelle nodes.
   mdofs = idofm(3)+[4:6];
   [Tmm0,dTmm0] = dTdth (q(idofm(3)+[4:6]));
   Tm0t = TFromTheta (P(idofm(3)+[4:6]));
   ir = 3*Ndj + 9*Nnod + mdofs;
   ic = iyy + 3;
   ddy(ir,ic) = -Tm0t*Tmm0(:,1);
   ic = mdofs;
   for ith = 1:3
      ic3 = 3*(ith-1);
      ddy(ir,ic(ith)) = ddy(ir,ic(ith)) ...
                      - Tm0t*dTmm0(:,ic3+1)*yyout;
   end
   ir = 3*Ndj + 9*Nnod + idofs(3) + 6;
   ic = iyy + 3;
   ddy(ir,ic) = 1;

   % Link the generator torque to the nacelle and driveshaft nodes.
   [Tmm0,dTmm0] = dTdth (qgen(10:12));
   Tm0y = TFromTheta (Pgen(10:12));
   [Tnn0,dTnn0] = dTdth (qgen(22:24));
   Tn0d = TFromTheta (Pgen(22:24));
   ir = 3*Ndj + 9*Nnod + igen(10:12);
   ic = iye + 2;
   ddy(ir,ic) = -Tm0y*Tmm0(:,1);
   ic = igen(10:12);
   for ith = 1:3
      ic3 = 3*(ith-1);
      ddy(ir,ic(ith)) = ddy(ir,ic(ith)) ...
                      - Tm0y*dTmm0(:,ic3+1)*yeout(1);
   end
   ir = 3*Ndj + 9*Nnod + igen(22:24);
   ic = iye + 2;
   ddy(ir,ic) = Tn0d*Tnn0(:,1);
   ic = igen(22:24);
   for ith = 1:3
      ic3 = 3*(ith-1);
      ddy(ir,ic(ith)) = ddy(ir,ic(ith)) ...
                      + Tn0d*dTnn0(:,ic3+1)*yeout(1);
   end

   % Link the generator rotor speed to the generator electrical speed input.
   ir = iye + 1;
   ic = [igen Ndj+igen];
   ddy(ir,ic) = 0.5*epar(1)*ddyw;

   % Implement the links and reduce the state space.
   %
   % Special logic is used to speed up the inverse of Dy.  We take advantage
   % of the fact that the aerodynamic variables are largely decoupled with
   % each other.  That is, the aerodynamic part is pretty close to diagonal,
   % and can be inverted quickly.  This operation is done separately from
   % the rest of the Dy matrix.
   adofs = 4*Ndj+9*Nnod+[1:137*Nae].';
   rdofs = [[1:4*Ndj+9*Nnod] [4*Ndj+9*Nnod+137*Nae+1:Ny]].';
   Na = size(adofs,1);
   Nr = size(rdofs,1);
   Nz = size(cc,2) + size(ddu,2);
   H = speye(Ny) - ddy;
   Hrr = H(rdofs,rdofs);
   Hra = H(rdofs,adofs);
   Har = H(adofs,rdofs);
   Haa = H(adofs,adofs);
   RR  = [cc ddu];
   Ra  = RR(adofs,:);
   Rr  = RR(rdofs,:);
   HHR = Haa\[Har Ra];
   H4  = Hrr - Hra*HHR(:,1:Nr);
   Sr  = H4\(Rr - Hra*HHR(:,Nr+[1:Nz]));
   Sa  = -HHR(:,1:Nr)*Sr + HHR(:,Nr+[1:Nz]);
   C  = spalloc(Ny,Nx,nnz(Sr(:,1:Nx))+nnz(Sa(:,1:Nx)));
   D  = spalloc(Ny,Nu,nnz(Sr(:,Nx+[1:Nu]))+nnz(Sa(:,Nx+[1:Nu])));
   C(rdofs,:) = Sr(:,1:Nx);
   C(adofs,:) = Sa(:,1:Nx);
   D(rdofs,:) = Sr(:,Nx+[1:Nu]);
   D(adofs,:) = Sa(:,Nx+[1:Nu]);
   A = aa + bby*C;
   B = bbu + bby*D;

else

   Nyp = size(ccp,1) - Nup;
   Nyy = size(ccy,1) - Nuy;
   Nye = size(yeout,1) + size(yein,1) - Nue;
   Ny = Nyal + Nyp + Nyy + Nye;

   A = sparse(Nx,Nx);
   B = sparse(Nx,Nu);
   C = sparse(Ny,Nx);
   D = sparse(Ny,Nu);

end

