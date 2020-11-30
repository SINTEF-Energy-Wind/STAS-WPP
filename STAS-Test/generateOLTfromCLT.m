
clear;

pkg load control;

[length,time,mass,current,voltage,        ...
 velocity,force,power,stress,ndens,nvisc, ...
 stiffness,damping,resistance,inductance, ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

nm = 'DTU10MW';
inpnm = '_P060_';
opnm = '_P060_';
convertToAeroModes = 0;

Vmag = 10;
thV = 0;
Vg0 = [Vmag*cos(thV); Vmag*sin(thV); 0];

yaw = 0;
betas = [0;0;0];
azi = 0;
cp0 = cos(azi);
sp0 = sin(azi);

Tpw = 6;
Hsw = 2;
wang = 0;
fw = 1/Tpw;
aw = 2*pi*fw;
zetw = 0.1;
atw = 0.01*(2*pi);
Fw0 = 1*[cos(wang), sin(wang)].';

aV = 0.05*(2*pi);
% a3P computed based on W.
zet3P = 0.1;

eval(['[s,a] = STASTurbine_'  nm ' ();']);
eval(['epar  = STASElectric_' nm ' ();']);
eval(['ppar  = STASPitch_'    nm ' ();']);
eval(['ypar  = STASYaw_'      nm ' ();']);
eval(['c     = STASControl_'  nm ' ();']);
eval(['m     = STASSensor_'   nm ' ();']);

a.dens = 1.225                     / ndens;
a.visc = 1.789e-5                  / nvisc;
grav   = [0;0;0];  % Consistent with asymWind, gravity is implemented as
                   % an external force rather than part of the matrices.

vs     = [33000;0]                 / voltage;    % Grid electrical voltage.
we     = 50*(2*pi)                 * time;       % Grid electrical frequency.
th_e   = 0;                                      % Ref. elec. angle.
Vhdc   = c.Vhdc;                                 % DC link voltage command.
Qh     = 0                         / power;      % Reactive power command.

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
[imdofs,Nmd] = getmdofRefs (s);
Ndj = Ndof + 6;
Nnod = Ndof/6;

Pin = assemblePin (s);
[q,P,Ts0_B,TB0_g] =                                                     ...
      undeformedPosition (Pin,yaw,s.nacelle.delta,azi,s.driveshaft.phi, ...
                          betas,0,idofs,idofm,inods,inodm);
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);

Nb  = a.Nb;
Neb = a.Neb;
Nae = Nb*Neb;

Nmud = s.foundation.Nmud;
Nwater = s.foundation.Nwater;

% Load the operating-point state vector and closed-loop matrices
% based on the default controller.
txt = inpnm;
if (Vmag >= 10)
   Vstr = ['V' int2str(10*Vmag)];
else
   Vstr = ['V0' int2str(10*Vmag)];
end

eval(["load 'xpsi" txt Vstr ".txt';"]);
eval(["load 'dxpsi" txt Vstr ".txt';"]);
eval(["load 'Rgrav" txt Vstr ".txt';"]);
eval(["load 'ypsi" txt Vstr ".txt';"]);
eval(["load 'upsi" txt Vstr ".txt';"]);
eval(["load 'Lpsi" txt Vstr ".bin';"]);
eval(["load 'Apsi" txt Vstr ".bin';"]);
eval(["load 'Bpsi" txt Vstr ".bin';"]);
eval(["load 'Cpsi" txt Vstr ".bin';"]);
eval(["load 'Dpsi" txt Vstr ".bin';"]);
eval(["load 'dret" txt Vstr ".bin';"]);
eval(["load 'shape" txt Vstr ".bin';"]);
eval(["load 'mdamp" txt Vstr ".bin';"]);
eval(["load 'bldof" txt Vstr ".bin';"]);
eval(["xpsi = xpsi" txt Vstr ";"]);
eval(["dxpsi = dxpsi" txt Vstr ";"]);
eval(["Rgrav = Rgrav" txt Vstr ";"]);
eval(["ypsi = ypsi" txt Vstr ";"]);
eval(["upsi = upsi" txt Vstr ";"]);

slv = slaveDOFs (idofs);
vec = [1:Ndj].';
[jnk,ret,jnk2] = partitionMatrix (vec,slv,[]);
Nret = size(ret,1);

Neta = size(shape0,2);
Nxs  = 2*Neta;
Nus  = Ndj + Neta;
Nys  = 4*Ndj;
Nxa  = 7*3*size(a.icp,1);
Nua  = 3*Nae;
Nya  = 0;
Nxb  = 6;
Nub  = 3;
Nyb  = 3;
Nxy  = 2;
Nuy  = 1;
Nyy  = 1;
Nxe  = 25;
Nue  = 8;
Nye  = 3;
Nxm  = 10;
Num  = 0;
Nym  = 10;
Nxc  = 3;
Nuc  = 1;
Nyc  = 2;
Nysl = 4*Ndj + 9*Nnod;
Nyal = 137*Nae + 5;
Nybl = 9;
Nyyl = 3;
Nyel = 4;

jxs  = 0;
jus  = 0;
jys  = 0;
jxa  = jxs + Nxs;
jua  = jus + Nus;
jya  = jys + Nys;
jxb  = jxa + Nxa;
jub  = jua + Nua;
jyb  = jya + Nya;
jxy  = jxb + Nxb;
juy  = jub + Nub;
jyy  = jyb + Nyb;
jxe  = jxy + Nxy;
jue  = juy + Nuy;
jye  = jyy + Nyy;
jxm  = jxe + Nxe;
jum  = jue + Nue;
jym  = jye + Nye;
jxc  = jxm + Nxm;
juc  = jum + Num;
jyc  = jym + Nym;

Nx = jxc + Nxc;
Nu = juc + Nuc;
Ny = jyc + Nyc;

iazi = Neta - 4;
iW   = 2*Neta - 4;
igen = [idofs(3)+[1:6] idofm(5)+[1:6] idofs(4)+[1:6] idofs(5)+[1:6]].';
ipit = [Ndof+[4:6]].';
iyaw = Ndof + 1;

% Retain structures, aero, actuators, and electrical systems
% from xpsi.  Add in sensor and gen P control states and build
% an appropriate input vector.
oret = [1:jxm].';
Pes = (xpsi(jxe+[17 18]).')*xpsi(jxe+[19 20]);
qB  = ypsi(idofs(2)+[1:6]);
qBd = ypsi(Ndj+idofs(2)+[1:6]);
PB  = P(idofs(2)+[1:6]);
qn  = ypsi(idofm(3)+[1:6]);
qnd = ypsi(Ndj+idofm(3)+[1:6]);
Pn  = P(idofm(3)+[1:6]);
[vnac,ddyv] = globalVelocity (qn,qB,Pn,PB,qnd,qBd);

xop  = zeros(Nx,1);
dxop = zeros(Nx,1);
uop  = zeros(Nu,1);

xgp = xpsi(jxm+[16 17 18]);   % From default controller.
xop(1:jxm)    = xpsi(oret);
xop(jxm+1:Nx) = [xpsi(2*Neta-4);       ...  % W
                 xpsi(Neta-[2 1 0]);   ...  % beta
                 xpsi(Neta-5);         ...  % yaw
                 Pes;                  ...
                 vnac(1:2);            ...
                 Vmag;                 ...
                 0;                    ...  % thw
                 xgp];

dxop(1:jxm) = dxpsi(oret);
dxop(jxm+1:Nx) = zeros(Nx-jxm,1);

uop = [upsi(1:Ndj); dxpsi(Neta+[1:Neta]); upsi(Ndj+[1:3*Nae]); ...
       xpsi(Neta-[2 1 0]); xpsi(Neta-5); we; th_e; vs;         ...
       xpsi(2*Neta+Nxa+8+1); 0; Vhdc; Qh; Pes];  % ighq will be determined by P control.

qpsi   = ypsi(1:Ndj);
qdpsi  = ypsi(Ndj+[1:Ndj]);
qddpsi = ypsi(2*Ndj+[1:Ndj]);
Fpsi   = ypsi(3*Ndj+[1:Ndj]);
[b1,b2,b3] = MBCindices_Ndj (Ndj,idofs);
[TpsiB_Ndj,TBpsi_Ndj] = MBC (Ndj,b1,b2,b3,azi);
q   = TpsiB_Ndj*qpsi;
qd  = TpsiB_Ndj*qdpsi;

% Aerodynamic parameters and modes.
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);
[Tas,ch,Lel,foilwt,aoaz,aoast,xas,yas,iq] = BEMsetup (s,a);

Td_n = [cp0 -sp0 0;sp0 cp0 0;0 0 1];
Try = Tn_y;
Tyy0 = TFromTheta (q(idofs(3)+[4:6]));
Ty0g = TFromTheta (P(idofs(3)+[4:6]));
Trg = Ty0g*Tyy0*Try;

[Tar,Trg,TB0g,TsB,TBB0,dTar,dTsB,dTBB0,wga] = ...
                      BEMprepTransforms (s,a,q,qd,P,Tas);

[zr,Area,Dp,rp,Lp,xeg,xhg,xyg] = ...
         BEMprepProjections (s,a,q,P,Try,Trg);

bsh = bladeModeShape (s,ret,shape0);
Psia = aeroPsi (a,rp,bsh);

% Conversion matrix from aero spline to mode.
if (convertToAeroModes == 1)

   abm = a;
   abm.icp = [-1, -3].';  % First two flap modes.
   Psiam = aeroPsi (abm,rp,bsh);
   Nxam = size(Psiam,2);
   Nxar = Nx - (Nxa-Nxam);

   % Psia needs to be in body coordinates, while to transform xop I need
   % MBC coordinates.  This will involve some transforming to and from
   % MBC, then.
   b1 = [1:Nxa/3].';
   b2 = Nxa/3 + [1:Nxa/3].';
   b3 = 2*Nxa/3 + [1:Nxa/3].';
   [TpsiBa,TBpsia] = MBC (Nxa,b1,b2,b3,azi);
   b1 = [1:Nxam/3].';
   b2 = Nxam/3 + [1:Nxam/3].';
   b3 = 2*Nxam/3 + [1:Nxam/3].';
   [TpsiBm,TBpsim] = MBC (Nxam,b1,b2,b3,azi);

   Tsm = sparse (Nxar,Nx);
   Tsm(1:jxa,1:jxa) = speye(Nxs);
   Tsm(jxa+[1:Nxam],jxa+1:jxb) = TBpsim*(Psiam.')*Psia*TpsiBa;
   Tsm(jxa+Nxam+[1:Nx-jxb],jxb+1:Nx) = speye(Nx-jxb);

   xop = Tsm*xop;
   dxop = Tsm*dxop;
   a = abm;
   Psia = Psiam;

   % Re-index.
   Nxa  = Nxam;
   Nx   = Nxar;
   jxb  = jxa + Nxa;
   jxy  = jxb + Nxb;
   jxe  = jxy + Nxy;
   jxm  = jxe + Nxe;
   jxc  = jxm + Nxm;

end

% Generate the open-loop turbine matrices in MBC coordinates.
% xop, dxop, uop are input in MBC.  Aero inputs including Psia
% are in body coordinates.
[Lop,Rop,yop,Aop,Bop,Cop,Dop] =                    ...
         MBCOLT (1,xop,dxop,uop,s,a,               ...
                 epar,ppar,ypar,m,c,               ...
                 grav,P,shape0,mdamp0,             ...
                 Tas,Try,ch,Lel,foilwt,aoaz,aoast, ...
                 xas,yas,Psia,igen,ipit,iyaw);

txt = opnm;

if (Vmag >= 10)
   Vstr = ['V' int2str(10*Vmag)];
else
   Vstr = ['V0' int2str(10*Vmag)];
end

eval(["save('-ascii', 'xop" txt Vstr ".txt','xop');"]);
eval(["save('-ascii', 'uop" txt Vstr ".txt','uop');"]);
eval(["save('-ascii','dxop" txt Vstr ".txt','dxop');"]);
eval(["save('-ascii', 'yop" txt Vstr ".txt','yop');"]);
eval(["save('-binary','Lop" txt Vstr ".bin','Lop');"]);
eval(["save('-binary','Aop" txt Vstr ".bin','Aop');"]);
eval(["save('-binary','Bop" txt Vstr ".bin','Bop');"]);
eval(["save('-binary','Cop" txt Vstr ".bin','Cop');"]);
eval(["save('-binary','Dop" txt Vstr ".bin','Dop');"]);

dret = [[1:Neta-6] Neta-[5 3 2 1 0] [(Neta+1):Nx]].'; 
[slap,shp,ifrq] = eigVal (Lop(dret,dret)\Aop(dret,dret));



