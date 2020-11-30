% Creates an observer model.
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
% ------------------------- Sensors ---------------------------
%   Wm         1         W          1
%   bm        2:4        bet       2:4
%   ym         5         yaw        5
%   Pem        6         Pe         6
%   vm        7,8        vnac      7,8
%   Vam        9         Va         9
%   thm       10         tha       10
%
% ------------------------ P control --------------------------
%  (Pem)                 ighq       1              Phat      1
%   PsiP       1         Pe         2             (Pe)
%   nd        2,3

clear;

pkg load control;

nm = 'DTU10MW_mod';
inpnm = '_mod_kf11_'; % '_mod_'; % '_asym_'; % '_mod_'; % '_idl_b14_W15_'; % 
opnm = '_mod_kf11_'; % '_mod_tf11_'; % '_'; % '_mod_tf11_'; % '_idl_b14_W15_gobs_'; % 
outnm = '_mod_kf11_'; % '_mod_tf11_'; % '_'; % '_mod_tf11_'; % '_idl_b14_W15_'; % 
genmat = 1;
xpsiflag = 1;
convertToAeroModes = 1;

Vmag = 10;
thV = 0;
Vg0 = [Vmag*cos(thV); Vmag*sin(thV); 0];
Tpw = 6;
Hsw = 2;
fw = 1/Tpw;
aw = 2*pi*fw;
zetw = 0.1;
atw = 0.01*(2*pi);
aV = 0.05*(2*pi);
% a3P computed based on W.
zet3P = 0.1;
wang = pi/2;
Fw0 = 1*[cos(wang), sin(wang)].';

yaw = 0;
betas = [0;0;0];
azi = 0;
cp0 = cos(azi);
sp0 = sin(azi);

Nmud = 5;

freqs = [[0.0001:0.0001:0.001] [0.002:0.001:0.1] [0.102:0.002:2] [2.1:0.1:10] [11:1:100]].'; % For plotting.
%freqs = [[0.005:0.005:4] [4.1:0.1:20] [21:100]].'; % For TF projections.
Nf = size(freqs,1);

eval(['[s,a] = STASTurbine_'  nm ' ();']);
eval(['epar  = STASElectric_' nm ' ();']);
eval(['ppar  = STASPitch_'    nm ' ();']);
eval(['ypar  = STASYaw_'      nm ' ();']);
eval(['c     = STASControl_'  nm ' ();']);
eval(['m     = STASSensor_'   nm ' ();']);

[length,time,mass,current,voltage,        ...
 velocity,force,power,stress,ndens,nvisc, ...
 stiffness,damping,resistance,inductance, ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

%=============================
%c.Kp    = -0.00015                / (current/power);  % -0.00015
%c.Ki    = -0.00250                 / (current/(power*time));  % -0.0025
%=============================


a.dens = 1.225                     / ndens;
a.visc = 1.789e-5                  / nvisc;
%grav   = [0;0;-9.807]              / (length/(time^2));
grav   = [0;0;0];  % Consider gravity constant, irrelevant for perturbed
                   % equations.

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
if (xpsiflag == 1)
   eval(["load 'xpsi" txt "V" int2str(Vmag) ".txt';"]);
   eval(["load 'dxpsi" txt "V" int2str(Vmag) ".txt';"]);
   eval(["load 'Rgrav" txt "V" int2str(Vmag) ".txt';"]);
   eval(["load 'ypsi" txt "V" int2str(Vmag) ".txt';"]);
   eval(["load 'upsi" txt "V" int2str(Vmag) ".txt';"]);
   eval(["load 'Lpsi" txt Vstr ".bin';"]);
   eval(["load 'Apsi" txt Vstr ".bin';"]);
   eval(["load 'Bpsi" txt Vstr ".bin';"]);
   eval(["load 'Cpsi" txt Vstr ".bin';"]);
   eval(["load 'Dpsi" txt Vstr ".bin';"]);
   eval(["load 'dret" txt Vstr ".bin';"]);
   eval(["load 'shape" txt Vstr ".bin';"]);
   eval(["load 'mdamp" txt Vstr ".bin';"]);
   eval(["load 'bldof" txt Vstr ".bin';"]);
   eval(["xpsi = xpsi" txt "V" int2str(Vmag) ";"]);
   eval(["dxpsi = dxpsi" txt "V" int2str(Vmag) ";"]);
   eval(["Rgrav = Rgrav" txt "V" int2str(Vmag) ";"]);
   eval(["ypsi = ypsi" txt "V" int2str(Vmag) ";"]);
   eval(["upsi = upsi" txt "V" int2str(Vmag) ";"]);
else
   eval(["load 'xop" txt "V" int2str(Vmag) ".txt';"]);
   eval(["load 'dxop" txt "V" int2str(Vmag) ".txt';"]);
   eval(["load 'Rgrav" txt "V" int2str(Vmag) ".txt';"]);
   eval(["load 'yop" txt "V" int2str(Vmag) ".txt';"]);
   eval(["load 'uop" txt "V" int2str(Vmag) ".txt';"]);
   eval(["load 'Lop" txt Vstr ".bin';"]);
   eval(["load 'Aop" txt Vstr ".bin';"]);
   eval(["load 'dret" txt Vstr ".bin';"]);
   eval(["load 'shape" txt Vstr ".bin';"]);
   eval(["load 'mdamp" txt Vstr ".bin';"]);
   eval(["load 'bldof" txt Vstr ".bin';"]);
   eval(["xop = xop" txt "V" int2str(Vmag) ";"]);
   eval(["dxop = dxop" txt "V" int2str(Vmag) ";"]);
   eval(["Rgrav = Rgrav" txt "V" int2str(Vmag) ";"]);
   eval(["yop = yop" txt "V" int2str(Vmag) ";"]);
   eval(["uop = uop" txt "V" int2str(Vmag) ";"]);
   xpsi = xop;
   dxpsi = dxop;
   ypsi = yop;
   upsi = uop;
end
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





%{
W = xop(iW);
wg = 0.5*W*epar(1);

xel = xop(jxe+[1:25]);
uel = [wg; uop(jue+[1:8])];
[dxeldt,yel,aae,bbe,cce,dde] = buildTurbineElectric (1,xel,uel,epar);

xgp = [xop(jxm+6); xop(jxc+[1:3])];
ugp = [uop(juc+1); uop(juc+1)];
gppar = [c.Kp c.Ki c.ap c.anp c.z1np c.z2np];
[dxgpdt,ygp,aag,bbg,ccg,ddg] = genPcontrol (xgp,ugp,gppar);

nnx = size(aae,1) + size(aag,1);
aa  = zeros(nnx,nnx);
bb  = zeros(nnx,1);
aa(1:Nxe,1:Nxe) = aae;
aa(Nxe+1:nnx,Nxe+1:nnx) = aag;
aa(1:Nxe,Nxe+1:nnx) = bbe(:,7)*ccg(1,:);
aa(Nxe+1:nnx,1:Nxe) = bbg(:,2)*vs(1)*cce(11,:);


%[slap,shp,ifrq] = eigVal (aae);
%[slap,shp,ifrq] = eigVal (aag);
[slap,shp,ifrq] = eigVal (aa);

%ist = [jxe+[1:25] jxm+6 jxc+[1:3]].';
%[slap,shp,ifrq] = eigVal (Lop(ist,ist)\Aop(ist,ist));

return

%}






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

if (xpsiflag == 1)

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

end

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



CAUTION, MAY BE AN ERROR ASSOCIATED WITH MBC OF AERO DOFS.  SEE GENERATEOLTFROMCLT.M



% Conversion matrix from aero spline to mode.
if (convertToAeroModes == 1)
   abm = a;
   abm.icp = [-1, -3].';  % First two flap modes.
   Psiam = aeroPsi (abm,rp,bsh);
   Nxam = size(Psiam,2);
   Nxar = Nx - (Nxa-Nxam);

   Tsm = sparse (Nxar,Nx);
   Tsm(1:jxa,1:jxa) = speye(Nxs);
   Tsm(jxa+[1:Nxam],jxa+1:jxb) = (Psiam.')*Psia;
   Tsm(jxa+Nxam+[1:Nx-jxb],jxb+1:Nx) = speye(Nx-jxb);

   xop = Tsm*xop;      % CAUTION, ISN'T PSIA IN BODY COORDINATES AND XOP IN MBC?
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

if (genmat == 1)

   % Generate the open-loop turbine matrices in MBC coordinates.
   [Lop,Rop,yop,Aop,Bop,Cop,Dop] =                    ...
            MBCOLT (1,xop,dxop,uop,s,a,               ...
                    epar,ppar,ypar,m,c,               ...
                    grav,P,shape0,mdamp0,             ...
                    Tas,Try,ch,Lel,foilwt,aoaz,aoast, ...
                    xas,yas,Psia,igen,ipit,iyaw);

   txt = opnm;
   eval(["save('-ascii','xop" txt "V" int2str(round(Vmag)) ".txt','xop');"]);
   eval(["save('-ascii','uop" txt "V" int2str(round(Vmag)) ".txt','uop');"]);
   eval(["save('-ascii','dxop" txt "V" int2str(round(Vmag)) ".txt','dxop');"]);
   eval(["save('-ascii','yop" txt "V" int2str(round(Vmag)) ".txt','yop');"]);

   if (Vmag < 10)
      eval(["save('-binary','Lop" txt "V0" int2str(round(10*Vmag)) ".bin','Lop');"]);
      eval(["save('-binary','Aop" txt "V0" int2str(round(10*Vmag)) ".bin','Aop');"]);
      eval(["save('-binary','Bop" txt "V0" int2str(round(10*Vmag)) ".bin','Bop');"]);
      eval(["save('-binary','Cop" txt "V0" int2str(round(10*Vmag)) ".bin','Cop');"]);
      eval(["save('-binary','Dop" txt "V0" int2str(round(10*Vmag)) ".bin','Dop');"]);
   else
      eval(["save('-binary','Lop" txt "V"  int2str(round(10*Vmag)) ".bin','Lop');"]);
      eval(["save('-binary','Aop" txt "V"  int2str(round(10*Vmag)) ".bin','Aop');"]);
      eval(["save('-binary','Bop" txt "V"  int2str(round(10*Vmag)) ".bin','Bop');"]);
      eval(["save('-binary','Cop" txt "V"  int2str(round(10*Vmag)) ".bin','Cop');"]);
      eval(["save('-binary','Dop" txt "V"  int2str(round(10*Vmag)) ".bin','Dop');"]);
   end

dret = [[1:Neta-6] Neta-[5 3 2 1 0] [(Neta+1):Nx]].'; 
[slap,shp,ifrq] = eigVal (Lop(dret,dret)\Aop(dret,dret));

return

else

   txt = opnm;
   eval(["load 'xop"   txt "V" int2str(Vmag) ".txt';"]);
   eval(["load 'dxop"  txt "V" int2str(Vmag) ".txt';"]);
   eval(["load 'yop"   txt "V" int2str(Vmag) ".txt';"]);
   eval(["load 'uop"   txt "V" int2str(Vmag) ".txt';"]);
   if (Vmag < 10)
      eval(["load 'Lop"   txt "V0" int2str(10*Vmag) ".bin';"]);
      eval(["load 'Aop"   txt "V0" int2str(10*Vmag) ".bin';"]);
      eval(["load 'Bop"   txt "V0" int2str(10*Vmag) ".bin';"]);
      eval(["load 'Cop"   txt "V0" int2str(10*Vmag) ".bin';"]);
      eval(["load 'Dop"   txt "V0" int2str(10*Vmag) ".bin';"]);
   else
      eval(["load 'Lop"   txt "V"  int2str(10*Vmag) ".bin';"]);
      eval(["load 'Aop"   txt "V"  int2str(10*Vmag) ".bin';"]);
      eval(["load 'Bop"   txt "V"  int2str(10*Vmag) ".bin';"]);
      eval(["load 'Cop"   txt "V"  int2str(10*Vmag) ".bin';"]);
      eval(["load 'Dop"   txt "V"  int2str(10*Vmag) ".bin';"]);
   end
   eval(["xop = xop"   txt "V" int2str(Vmag) ";"]);
   eval(["dxop = dxop" txt "V" int2str(Vmag) ";"]);
   eval(["yop = yop"   txt "V" int2str(Vmag) ";"]);
   eval(["uop = uop"   txt "V" int2str(Vmag) ";"]);

end

% Eliminate the azimuth DOF.
dret = [[1:Neta-6] Neta-[5 3 2 1 0] [(Neta+1):Nx]].'; 
%dret = [[1:Neta-6] Neta+[1:Neta-6] [(2*Neta+1):Nx]].'; 
Ndr = size(dret,1);

% Get the open-loop eigenmodes.  Always partition/eliminate first,
% then invert L afterwards.
LA = Lop(dret,dret)\Aop(dret,dret);
%[slap,shp,ifrq] = eigVal (LA);
[slap,shp,ifrq] = eigVal_silent (LA);
parfac;

% Diagonalize.
Phim = shp(:,ifrq);
Psim = inv(Phim);
ii = [1:Ndr].';
jj = ii;
Lam = sparse (ii,jj,slap,Ndr,Ndr);

% Some of the very closely spaced eigenfrequencies may lie out of order.
% This causes problems when I assume that conjugate modes are mirrored
% in the sorted mode shape matrix.  The solution is to force complex
% conjugacy by copying half the complex mode shapes and eigenvalues to
% the other half.  This is done in the 'for' loop below.
iY = speye(Ndr);
Nosc = 0;
for im1 = 1:floor(Ndr/2)

   imn = Ndr - (im1-1);

   if (abs(imag(slap(im1))) > 0)
      Nosc = Nosc + 1;
      Phim(:,imn) = conj(Phim(:,im1));
      Psim(imn,:) = conj(Psim(im1,:));
      slap(imn) = conj(slap(im1));
      iY([im1 imn],[im1 imn]) = [1 1;-i i];
   end

end
Nexp = Ndr - 2*Nosc;
Nmds = Nosc + Nexp;

%=========================================================================
% Open-loop control-to-sensor transfer functions; these are the TFs that
% must be captured with the observer.  Also the disturbance-to-sensor TFs,
% if I want to observe the wind and waves.
% 
% Controls: b0, bc, bs, Phat, yawd
% Disturbances: Veff, thV, Vc, Vs, Fw, thw
% Sensors: W, beta, yaw, Pe, vnac^g(x,y), Va, tha^g
%
% There are 11 x 10 = 110 total TFs.
%=========================================================================

iuc = [jub+[1 2 3], juc+1, jub+4];
iwx0 = jua + [1:3:Nae-2].';
iwy0 = jua + [2:3:Nae-1].';
iwz0 = jua + [3:3:Nae].';
iwxc = jua + Nae + [1:3:Nae-2].';
iwyc = jua + Nae + [2:3:Nae-1].';
iwzc = jua + Nae + [3:3:Nae].';
iwxs = jua + 2*Nae + [1:3:Nae-2].';
iwys = jua + 2*Nae + [2:3:Nae-1].';
iwzs = jua + 2*Nae + [3:3:Nae].';
iFwx = 6*(Nmud+Nwater-1) + 1;
iFwy = 6*(Nmud+Nwater-1) + 2;

ixm = jxm + [1:Nxm].';
Nsens = size(ixm,1);

PLB0  = Psim*(Lop(dret,dret)\Bop(dret,:));
PLBu = PLB0(:,iuc);

% Create B matrices for the wind and wave disturbances in terms of
% amplitudes and directions.  These B matrices should accept as input
% components in (mag,ph) and return (x,y).
%
% Order is Vamp, thV, Vc, Vs, Fwamp, thw.
PLBw = sparse(Ndr,6);
mag = sqrt(Vg0(1)^2 + Vg0(2)^2);
ang = atan2(Vg0(2),Vg0(1));
ca = cos(ang);
sa = sin(ang);
PLBw(:,1:2) = [sum(PLB0(:,iwx0),2) sum(PLB0(:,iwy0),2)]*[ca, -mag*sa;sa, mag*ca];
PLBw(:,3:4) = [sum(PLB0(:,iwxc),2) sum(PLB0(:,iwxs),2)];
mag = sqrt(Fw0(1)^2 + Fw0(2)^2);
ang = atan2(Fw0(2),Fw0(1));
ca = cos(ang);
sa = sin(ang);
PLBw(:,5:6) = [PLB0(:,iFwx) PLB0(:,iFwy)]*[ca, -mag*sa;sa, mag*ca];

PLB = [PLBu PLBw];

CC = sparse(Nsens,Nx);
CC(:,ixm) = speye(Nsens);
CP = CC(:,dret)*Phim;

Ninp = 11;
Nout = Nsens;
NTFs = Ninp*Nout;
Hmod = zeros(Nf,NTFs*Ndr);

for ifreq = 1:Nf

%if (mod(ifreq,100) == 0)
%printf('%6d of %6d\n',ifreq,Nf);
%fflush(stdout);
%end

   f = freqs(ifreq);
   w = 2*pi*f;

   % Disturbances-to-states.
   iwslap = i*w - slap.';

   for jinp = 1:Ninp
      for jout = 1:Nout
         icn = Ndr*(jout-1) + Ndr*Nout*(jinp-1);
         Hmod(ifreq,icn+[1:Ndr]) = CP(jout,:).*(PLB(:,jinp).')./iwslap;
      end
   end

end

H = zeros (Nf,NTFs);
for itf = 1:NTFs
   icn = Ndr*(itf-1);
   H(:,itf) = sum(Hmod(:,icn+[1:Ndr]),2);
end

%eval(["save ('-binary','Hmod_V" int2str(10*Vmag) ".bin','Hmod');"]);
%eval(["save ('-binary','H_V" int2str(10*Vmag) ".bin','H');"]);

Hmag = abs(H);
Hang = atan2(imag(H),real(H))/pi;

% Compute projections.  Use the max error in the projection up to 1 Hz
% as a metric.
proj = zeros(Nf,NTFs*Nmds);
for itf = 1:NTFs

   icr = Ndr*(itf-1);
   icn = Nmds*(itf-1);
   hvec = H(:,itf)./(conj(H(:,itf)).*H(:,itf));

   proj(:,icn+[1:Nosc]) = (conj(Hmod(:,icr+[1:Nosc]))             ...
                        +  conj(Hmod(:,icr+[Ndr:-1:Ndr-Nosc+1]))).*hvec;
   proj(:,icn+Nosc+[1:Nexp]) = conj(Hmod(:,icr+Nosc+[1:Nexp])).*hvec;

end

pval = zeros(Nf,NTFs*Nmds);
sindx = zeros(Nf,NTFs*Nmds);
for itf = 1:NTFs
   icn = Nmds*(itf-1);
   [pval(:,icn+[1:Nmds]),sindx(:,icn+[1:Nmds])] = ...
                       sort (abs(real(proj(:,icn+[1:Nmds]))),2,'descend');
   for ifreq = 1:Nf
      pval(ifreq,icn+[1:Nmds]) = proj(ifreq,icn+sindx(ifreq,icn+[1:Nmds]));
   end
end

pcum = zeros(Nf,NTFs*Nmds);
for itf = 1:NTFs

   icn = Nmds*(itf-1);

   pcum(:,icn+1) = pval(:,icn+1);
   for imd = 2:Nmds
      pcum(:,icn+imd) = pcum(:,icn+imd-1) + pval(:,icn+imd);
   end

end

writeOLTFTable;

% Reduced criteria for choosing modes.
%jf = [4 34 48 80].';  % 0.02, 0.17, 0.24, 0.40 Hz.
jf = [29 144 179 259].'; % ... same, for the higher-res frequencies for plotting.

%{
tret = [1 1 0 0 1 1 1 1 0 0  ...
        1 0 1 1 1 1 1 1 0 0  ...
        1 0 1 1 1 1 1 1 0 0  ...
        1 0 0 0 0 1 1 1 0 0  ...
        1 1 1 1 1 1 1 1 0 0  ...
        1 0 0 0 0 1 1 1 1 0  ...
        1 0 0 0 0 1 1 1 0 1  ...
        1 0 0 0 1 0 0 0 0 0  ...
        1 0 0 0 1 0 1 1 0 0  ...
        1 0 0 0 0 1 1 1 0 0  ...
        1 0 0 0 0 1 1 1 0 0; ...
%       ---------
        1 1 0 0 1 1 1 1 0 0  ...
        1 0 1 1 1 1 1 1 0 0  ...
        1 0 1 1 1 1 1 1 0 0  ...
        1 0 0 0 0 1 1 1 0 0  ...
        1 1 1 1 1 1 1 1 0 0  ...
        1 0 0 0 0 1 1 1 1 0  ...
        1 0 0 0 0 1 1 1 0 1  ...
        0 0 0 0 1 0 1 0 0 0  ...
        1 0 0 0 0 0 1 1 0 0  ...
        1 0 0 0 0 1 1 1 0 0  ...
        1 0 0 0 0 1 1 1 0 0; ...
%       ---------
        1 1 0 0 1 1 1 1 0 0  ...
        1 0 1 1 1 1 1 1 0 0  ...
        1 0 1 1 1 1 1 1 0 0  ...
        1 0 0 0 0 1 1 1 0 0  ...
        1 1 1 1 1 1 1 1 0 0  ...
        1 0 0 0 0 1 1 1 1 0  ...
        1 0 0 0 0 1 1 1 0 1  ...
        1 0 0 0 0 1 1 1 0 0  ...
        1 0 0 0 1 1 1 1 0 0  ...
        1 0 0 0 0 1 1 1 0 0  ...
        1 0 0 0 0 1 1 1 0 0; ...
%       ---------
        1 1 0 0 1 1 1 1 0 0  ...
        1 0 1 1 1 1 1 1 0 0  ...
        1 0 1 1 1 1 1 1 0 0  ...
        1 0 0 0 0 1 1 1 0 0  ...
        1 1 1 1 1 1 1 1 0 0  ...
        1 0 0 0 0 1 1 1 1 0  ...
        1 0 0 0 0 1 1 1 0 1  ...
        1 0 0 0 0 1 1 1 0 0  ...
        1 0 0 0 1 1 1 1 0 0  ...
        1 0 0 0 0 1 1 1 0 0  ...
        1 0 0 0 0 1 1 1 0 0];
%}

tret = [1 1 1 1 1 0 1 0 0 0  ...
        1 1 1 1 0 0 1 1 0 0  ...
        1 1 1 1 0 0 1 1 0 0  ...
        1 0 0 0 1 1 1 1 0 0  ...
        1 0 0 0 1 0 1 1 0 0  ...
        1 0 0 0 0 0 1 1 1 0  ...
        1 0 0 0 0 0 0 0 0 1  ...
        1 0 0 0 0 0 0 0 0 0  ...
        1 0 0 0 0 0 1 1 0 0  ...
        0 0 0 0 0 0 0 1 0 0  ...
        0 0 0 0 0 0 1 0 0 0; ...
%       ---------
        1 1 1 1 1 0 1 0 0 0  ...
        1 1 1 1 0 0 1 1 0 0  ...
        1 1 1 1 0 0 1 1 0 0  ...
        1 0 0 0 1 1 1 1 0 0  ...
        1 0 0 0 1 0 1 1 0 0  ...
        1 0 0 0 0 0 1 1 1 0  ...
        1 0 0 0 0 0 0 0 0 1  ...
        1 0 0 0 0 0 0 0 0 0  ...
        1 0 0 0 0 0 1 1 0 0  ...
        0 0 0 0 0 0 0 1 0 0  ...
        0 0 0 0 0 0 1 0 0 0; ...
%       ---------
        1 1 1 1 1 0 1 0 0 0  ...
        1 1 1 1 0 0 1 1 0 0  ...
        1 1 1 1 0 0 1 1 0 0  ...
        1 0 0 0 1 1 1 1 0 0  ...
        1 0 0 0 1 0 1 1 0 0  ...
        1 0 0 0 0 0 1 1 1 0  ...
        1 0 0 0 0 0 0 0 0 1  ...
        1 0 0 0 0 0 0 0 0 0  ...
        1 0 0 0 0 0 1 1 0 0  ...
        0 0 0 0 0 0 0 1 0 0  ...
        0 0 0 0 0 0 1 0 0 0; ...
%       ---------
        1 1 1 1 1 0 1 0 0 0  ...
        1 1 1 1 0 0 1 1 0 0  ...
        1 1 1 1 0 0 1 1 0 0  ...
        1 0 0 0 1 1 1 1 0 0  ...
        1 0 0 0 1 0 1 1 0 0  ...
        1 0 0 0 0 0 1 1 1 0  ...
        1 0 0 0 0 0 0 0 0 1  ...
        1 0 0 0 0 0 0 0 0 0  ...
        1 0 0 0 0 0 1 1 0 0  ...
        0 0 0 0 0 0 0 1 0 0  ...
        0 0 0 0 0 0 1 0 0 0];

% Get a list of all the modes that together satisfy the above TFs
% to sufficient accuracy.
acc = 0.05;

NL = sum(tret(1,:));
pcrL = zeros(Nmds,NL);
sirL = zeros(Nmds,NL);
jj = 0;
for itf = 1:NTFs
   icn = Nmds*(itf-1);
   if (tret(1,itf) == 1)
      jj = jj + 1;
      pcrL(:,jj) = pcum(jf(1),icn+[1:Nmds]).';
      sirL(:,jj) = sindx(jf(1),icn+[1:Nmds]).';
   end
end

NW = sum(tret(2,:));
pcrW = zeros(Nmds,NW);
sirW = zeros(Nmds,NW);
jj = 0;
for itf = 1:NTFs
   icn = Nmds*(itf-1);
   if (tret(2,itf) == 1)
      jj = jj + 1;
      pcrW(:,jj) = pcum(jf(2),icn+[1:Nmds]).';
      sirW(:,jj) = sindx(jf(2),icn+[1:Nmds]).';
   end
end

NR = sum(tret(3,:));
pcrR = zeros(Nmds,NR);
sirR = zeros(Nmds,NR);
jj = 0;
for itf = 1:NTFs
   icn = Nmds*(itf-1);
   if (tret(3,itf) == 1)
      jj = jj + 1;
      pcrR(:,jj) = pcum(jf(3),icn+[1:Nmds]).';
      sirR(:,jj) = sindx(jf(3),icn+[1:Nmds]).';
   end
end

N3 = sum(tret(4,:));
pcr3 = zeros(Nmds,N3);
sir3 = zeros(Nmds,N3);
jj = 0;
for itf = 1:NTFs
   icn = Nmds*(itf-1);
   if (tret(4,itf) == 1)
      jj = jj + 1;
      pcr3(:,jj) = pcum(jf(4),icn+[1:Nmds]).';
      sir3(:,jj) = sindx(jf(4),icn+[1:Nmds]).';
   end
end

pcr = [pcrL pcrW pcrR pcr3];
sir = [sirL sirW sirR sir3];

NLWR = NL + NW + NR + N3;
mlist = zeros(Nmds*NLWR,1);
jj = 0;
for ii = 1:NLWR

   flg = 1;
   kk = 1;
   jj = jj + 1;
   mlist(jj) = sir(kk,ii);           % Always keep the first mode.
   e1 = abs(1 - pcr(kk,ii));
   e2 = 1;
   while (flg)
      kk = kk + 1;
      e2 = abs(1 - pcr(kk,ii));
      if ((e1 < acc) && (e2 < acc))
         flg = 0;
      else
         jj = jj + 1;
         mlist(jj) = sir(kk,ii);
         e1 = e2;
      end
   end

end

[mls, inds] = sort(mlist);
mls = mls(mls > 0);
Nmls = size(mls,1);

ml = zeros(Nmls,1);
kk = 0;
prev = 0;
for jj = 1:Nmls
   if (mls(jj) ~= prev)
      kk = kk + 1;
      ml(kk) = mls(jj);
      prev = ml(kk);
   end
end

ml = ml(ml > 0);
Nml = size(ml,1);

% Get the full set of indices (+ and - frequencies).
mlf = zeros(2*Nml,1);
kk = 0;
for jj = 1:Nml
   if (ml(jj) <= Nosc)
      kk = kk + 1;
      mlf(kk) = ml(jj);
      kk = kk + 1;
      mlf(kk) = Ndr - ml(jj) + 1;
   else
      kk = kk + 1;
      mlf(kk) = ml(jj);
   end
end

[mlf,indm] = sort (mlf);
mlf = mlf(mlf > 0);
Nxi = size(mlf,1);

% Build the initial observer matrices.
Nxaug = 7;     % 5: Vmag, thV, s3, V3, sw, Fw, thw.
Nxob = Nxi + Nxaug;
Nydis = 5;
Aoi = sparse(Nxob,Nxob);
Boi = sparse(Nxob,5);
Coi = sparse(Nsens+Nydis,Nxob);  % Include observed Vmag, thV, V3, Fw, thw
                                 % such that these are available after DOF reduction.

mfrq = slap(mlf);
Aoi(1:Nxi,1:Nxi) = sparse(diag(mfrq));
Boi(1:Nxi,:) = PLBu(mlf,:);             % b0,c,s, Phat, yaw.
Aoi(1:Nxi,Nxi+[1 2 4 6 7]) = PLBw(mlf,[1 2 1 5 6]); % V,thV,[s3],V3,[sw],Fw,thw.
Coi(1:Nsens,1:Nxi) = CP(:,mlf);         % Wm, bet, yaw, Pe, vnacx,y, Va, tha

% Augment with the environmental dynamics.  Wind: LP filter, 3P BP filter.
% Wave amplitude: sinusoidal, band-pass.  Wave direction: LP filter.
Aoi(Nxi+[1:2],Nxi+[1:2]) = [-aV, 0; ...
                             0, -aV];

W = xpsi(iW);
a3P = 3*W;
a32 = a3P^2;
tza3 = 2*zet3P*a3P;
Aoi(Nxi+[3:4],Nxi+[3:4]) = [0,     1; ...
                           -a32, -tza3];

aw2 = aw^2;
tza = 2*zetw*aw;
Aoi(Nxi+[5:6],Nxi+[5:6]) = [0,     1; ...
                           -aw2, -tza];
Aoi(Nxi+7,Nxi+7) = -atw;

% Add the observed states to the outputs.  These are not a part of the
% observer gain calculations.  They must be removed from the C matrix
% before computing the gains via LQR.  But they are added here so that
% they can be included as outputs to the controller, from the reduced-
% order observer model, where the states may be transformed beyond
% recognition.
Coi(Nsens+[1:Nydis],Nxi+[1 2 4 6 7]) = speye(5);  

% Perform the Y transform to get real matrices.
iYY = sparse(Nxob,Nxob);
iYY(1:Nxi,1:Nxi) = iY(mlf,mlf);
iYY(Nxi+[1:Nxaug],Nxi+[1:Nxaug]) = speye(Nxaug);
YY = inv(iYY);
Ayi = real(iYY*Aoi*YY);
Byi = real(iYY*Boi);
Cyi = real(Coi*YY);

Bwi = sparse(Nxob,6);  % Matrix to apply disturbances to the system model.
Bwi(1:Nxi,:) = real(iYY(1:Nxi,1:Nxi)*PLBw(mlf,:));

Bdi = sparse(Nxob,5);  % Matrix to inject noise to the environmental models.
Bdi(Nxi+[1:Nxaug],:) = [aV  0   0   0   0; ...  % Vmag
                        0  aV   0   0   0; ...  % thV
                        0   0   0   0   0; ...  % s3P
                        0   0  tza3 0   0; ...  % V3P
                        0   0   0   0   0; ...  % sw
                        0   0   0  tza  0; ...  % Fw
                        0   0   0   0  atw];    % thw

%-------------------
%dret = [1:Nxob-1].';   % Cut out the wave direction observation.
%Ayi = Ayi(dret,dret);
%Byi = Byi(dret,:);
%Bwi = Bwi(dret,:);
%Cyi = Cyi(:,dret);
%Nxob = size(dret,1);
%Nxaug = 4;
%-------------------

Aobs = Ayi;
Bobs = Byi;
Bdis = [Bwi Bdi];
Cobs = Cyi(1:Nsens,:);                  % Anemometer.
%Cobs = [Cyi(1:Nsens-2,:); zeros(2,Nxob)];  % Observation.
%Cobs(Nsens-1,Nxi+1) = 1;
%Cobs(Nsens,Nxi+2) = 1;
Cdis = Cyi(Nsens+[1:Nydis],:);
Dobs = sparse(size(Cobs,1),9);
Ddis = sparse(size(Cdis,1),9);

ret = [1:Nxob].';

Nxr = size(Aobs,1);

% Q:      Confidence in disturbance.
% R^(-1): Confidence in measurement.

Nmd = Nxr - Nxaug;

% Full-mode tuning.
obsflg = 1;
if (obsflg == 1)
   alf = 1;
   bet = 0.025;
   wR = bet*[1.0;   ... % W 
             0.052;   ... % b0
             0.052;   ... % bc
             0.052;   ... % bs
             0.2;   ... % yaw
             0.52;   ... % Pe
             0.52;   ... % vnx
             0.52;   ... % vny
            20.0;   ... % Va
             2.0];      % tha
   %    |--------controls---------||-sys disturbances-||--env states--|
%  wQ = [0.05 0.05 0.05 0.50 0.05    1 0.1 1 1 1 0.1    1 0.1 1 1 0.1].';
%  wQ = [0.05 0.05 0.05 0.50 0.05                       1 0.1 1 1 0.1].';
   wQ = [                            0  0  1 1 0  0     1 0.1 1 1 0.1].';
%  wQ = [                                               1 0.1 1 1 0.1].';
%  wQ = [                            1 0.1 1 1 1 0.1                 ].';
%   fqs = [0.02 0.17 0.24].';
%   wts = [0.33 0.33 0.33].';
   fqs = [0.02 0.17 0.24 0.40].';
   wts = [0.49 0.17 0.17 0.17].';

   Nfqs = size(fqs,1);
%   Bcomb = [Bobs Bdis(:,1:6)];    % Inject noise to controls + system.
%   Bcomb = [Bobs Bdis(:,7:11)];    % Inject noise to controls + env. model states.
%   Bcomb = [Bobs Bdis]; 
   Bcomb = Bdis;                    % Inject noise to system + env. model states.
%   Bcomb = Bdis(:,7:11);           % Inject noise to env. model states.
%   Bcomb = Bdis(:,1:6);            % Inject noise to system.
   TF = zeros (Nxr,size(Bcomb,2));
   for ifq = 1:Nfqs
      TF = TF + ((i*2*pi*fqs(ifq)*speye(Nxr) - Aobs)\Bcomb)*wts(ifq);
   end
   MM = TF;
   QQ = real(alf*MM*diag(wQ.^2)*(MM'));

end

RR = diag(wR.^2);

% There are zero eigenmodes, add a tiny stability margin.
%Aobs = Aobs - 1e-6*speye(Nxr);
%Aobs(Nxr,Nxr) = Aobs(Nxr,Nxr) - 1e-6;
%[slap,shp,ifrq] = eigVal(Aobs);

tf = 0.1;
dt = 0.0001;
KK0 = 0;
[PP,KK] = solveRiccati (1,2,Aobs,Cobs,QQ,RR,KK0,tf,dt);

KK0 = KK;
%KK0 = zeros(size(Cobs,1),size(Aobs,1));
mat = Aobs - (KK0.')*Cobs;
%[slap,shp,ifrq] = eigVal(mat);

[PP,KK] = solveRiccati (3,2,Aobs,Cobs,QQ,RR,KK0,tf,dt);

GG = KK.';
Aof = Aobs - GG*Cobs;
[slap,shp,ifrq] = eigVal_silent(Aof);
max(real(slap))

TFyy = zeros(Nf,100);  % 10 ystar's in and out.
TFdy = zeros(Nf,50);   % 10 ystar's in, 5 observed disturbances out.
TFyu = zeros(Nf,50);   % 5 controls in, 10 ystar's out.
TFdu = zeros(Nf,25);   % 5 controls in, 5 observed disturbances out.
for ifreq = 1:Nf
   w = 2*pi*freqs(ifreq);
   dxdy = (i*w*speye(Nxr) - Aof)\GG;
   dxdu = (i*w*speye(Nxr) - Aof)\Bobs;
   TFyy(ifreq,:) = reshape(Cobs*dxdy,100,1);
   TFdy(ifreq,:) = reshape(Cdis*dxdy,50,1);
   TFyu(ifreq,:) = reshape(Cobs*dxdu,50,1);
   TFdu(ifreq,:) = reshape(Cdis*dxdu,25,1);
end

figure(1);
clf;
set (gcf,'papersize',[8. 11.],'paperorientation','portrait', ...
     'paperposition',[0.1 0.2 8. 11.]);
clrz = [0    0    0; ...
        1    0    0; ...
        0.65 0    0; ...
        0.35 0    0; ...
        0    1    0; ...
        0    0.4  0; ...
        0    0    1; ...
        0.3  0.3  0.6; ...
        0.4  0.4  0.4; ...
        0.7  0.7  0.7];
posz = [0.05 0.825 0.42 0.15; ...
        0.05 0.630 0.42 0.15; ...
        0.05 0.435 0.42 0.15; ...
        0.05 0.240 0.42 0.15; ...
        0.05 0.045 0.42 0.15; ...
        0.55 0.825 0.42 0.15; ...
        0.55 0.630 0.42 0.15; ...
        0.55 0.435 0.42 0.15; ...
        0.55 0.240 0.42 0.15; ...
        0.55 0.045 0.42 0.15];

ylimz = [-6 1; -6 1; -6 1; -6 1; -6 1; -6 1; -6 1; -6 1; -6 1; -6 1; -6 1];
ypozz = [1 1 1 1 1 1 1 1 1 1];
celtxt = {'\Omega','\it\beta_{0}','\it\beta_{c}','\it\beta_{s}','\chi', ...
          '\itP_{e}','(\bfv_{\itn}\rm)_{\itX}','(\bfv_{\itn}\rm)_{\itY}','\itV_{\infty}', ...
          '\it\theta_{V}'};
hax = zeros(10,1);
for iax = 1:10
   hax(iax) = axes();
   set(hax(iax),'position',posz(iax,:),'xlim',[-4 2],'ylim',ylimz(iax,:), ...
       'fontsize',11,'fontname','Times-Roman','ticklength',[0 0], ...
       'xtick',[-4:1:2],'linewidth',1);
   xl = get(hax(iax),'xlim');
   yl = get(hax(iax),'ylim');
   text(xl(1)+1,yl(1)+ypozz(iax),celtxt{iax},'fontsize',11,'fontname','Times-Roman');
end

for iax = 1:10

   ic10 = 10*(iax-1);

   axes(hax(iax));
   hold on;
   box on;
   for jj = 1:10
      plot (log10(freqs),log10(abs(TFyy(:,ic10+jj))), ...
            'linewidth',1.2,'color',clrz(jj,:));
   end
   
end

saveas (gcf,'TFyy.pdf','pdf');

figure(2);
clf;
set (gcf,'papersize',[6. 8.],'paperorientation','landscape', ...
     'paperposition',[0.1 0.2 8. 6.]);
clrz = [0    0    0; ...
        1    0    0; ...
        0.65 0    0; ...
        0.35 0    0; ...
        0    1    0; ...
        0    0.4  0; ...
        0    0    1; ...
        0.3  0.3  0.6; ...
        0.4  0.4  0.4; ...
        0.7  0.7  0.7];
posz = [0.05 0.684 0.4 0.267; ...
        0.05 0.367 0.4 0.267; ...
        0.05 0.050 0.4 0.267; ...
        0.55 0.684 0.4 0.267; ...
        0.55 0.367 0.4 0.267];
ylimz = [-6 1; -6 1; -6 1; -6 1; -6 1];
ypozz = [1 1 1 1 1];
celtxt = {'\it\beta_{0}','\it\beta_{c}','\it\beta_{s}','\itP_{e}','\chi'};
hax = zeros(5,1);
for iax = 1:5
   hax(iax) = axes();
   set(hax(iax),'position',posz(iax,:),'xlim',[-4 2],'ylim',ylimz(iax,:), ...
       'fontsize',11,'fontname','Times-Roman','ticklength',[0 0], ...
       'xtick',[-4:1:2],'linewidth',1);
   xl = get(hax(iax),'xlim');
   yl = get(hax(iax),'ylim');
   text(xl(1)+1,yl(1)+ypozz(iax),celtxt{iax},'fontsize',11,'fontname','Times-Roman');
end

for iax = 1:5

   ic10 = 10*(iax-1);

   axes(hax(iax));
   hold on;
   box on;
   for jj = 1:10
      plot (log10(freqs),log10(abs(TFyu(:,ic10+jj))), ...
            'linewidth',1.2,'color',clrz(jj,:));
   end
   
end

saveas (gcf,'TFyu.pdf','pdf');

%{

figure(3);
clf;
for im = 1:10

   ic5  = 5*(im-1);

   subplot(2,5,im);
   loglog(freqs,abs(TFdy(:,ic5+1)), 'color',[0   0   0], ...
          freqs,abs(TFdy(:,ic5+2)), 'color',[0.5   0.5 0.5], ...
          freqs,abs(TFdy(:,ic5+3)), 'color',[0.4 0 0.4], ...
          freqs,abs(TFdy(:,ic5+4)), 'color',[0   0   1], ...
          freqs,abs(TFdy(:,ic5+5)), 'color',[0 0   0.5]);

end

%}
%{

figure(4);
clf;
for im = 1:5

   ic5  = 5*(im-1);

   subplot(2,3,im);
   loglog(freqs,abs(TFdy(:,ic5+1)), 'color',[0   0   0], ...
          freqs,abs(TFdy(:,ic5+2)), 'color',[0.5 0.5 0.5], ...
          freqs,abs(TFdy(:,ic5+3)), 'color',[0.4 0 0.4], ...
          freqs,abs(TFdy(:,ic5+4)), 'color',[0   0   1], ...
          freqs,abs(TFdy(:,ic5+5)), 'color',[0   0   0.5]);

end

%}

txt = outnm;
if (Vmag < 10)
   eval(["save('-binary','GG"   txt "V0" int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','GG');"]);
   eval(["save('-binary','Aobs" txt "V0" int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','Aobs');"]);
   eval(["save('-binary','Bobs" txt "V0" int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','Bobs');"]);
   eval(["save('-binary','Bdis" txt "V0" int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','Bdis');"]);
   eval(["save('-binary','Cobs" txt "V0" int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','Cobs');"]);
   eval(["save('-binary','Cdis" txt "V0" int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','Cdis');"]);
   eval(["save('-binary','Dobs" txt "V0" int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','Dobs');"]);
   eval(["save('-binary','Ddis" txt "V0" int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','Ddis');"]);
else
   eval(["save('-binary','GG"   txt "V"  int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','GG');"]);
   eval(["save('-binary','Aobs" txt "V"  int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','Aobs');"]);
   eval(["save('-binary','Bobs" txt "V"  int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','Bobs');"]);
   eval(["save('-binary','Bdis" txt "V"  int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','Bdis');"]);
   eval(["save('-binary','Cobs" txt "V"  int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','Cobs');"]);
   eval(["save('-binary','Cdis" txt "V"  int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','Cdis');"]);
   eval(["save('-binary','Dobs" txt "V"  int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','Dobs');"]);
   eval(["save('-binary','Ddis" txt "V"  int2str(round(10*Vmag)) "_wang" int2str(round((180/pi)*wang)) ".bin','Ddis');"]);
end



