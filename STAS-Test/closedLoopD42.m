
clear;

pkg load control;

[length,time,mass,current,voltage,              ...
 velocity,force,power,stress,density,viscosity, ...
 stiffness,damping,resistance,inductance,       ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

nm = 'DTU10MW';
inpnm = '_P060_';
convertToAeroModes = 0;

Vmag = 10                                   / velocity;
thV = 0;
Vg0 = [Vmag*cos(thV); Vmag*sin(thV); 0];
dens = 1.225                                / density;
Area = pi*(89.15^2)                         / length^2;
CPstar = 0.48;

yaw = 0;
betas = [0;0;0];
azi = 0;
cp0 = cos(azi);
sp0 = sin(azi);

Tpw = 6                                     / time;
Hsw = 2                                     / length;
wang = 0;
fw = 1/Tpw;
aw = 2*pi*fw;
zetw = 0.1;
atw = 0.01*(2*pi)                           * time;
Fw0 = 1.e6*[cos(wang), sin(wang)].'         / force;

aV = 0.05*(2*pi);
% a3P computed based on W.
zet3P = 0.1;

fco = 2.0                                   * time;  % > fco solved as static.
dt = 0.05                                   / time;  % Set to -1 for auto.

Nf = 2^12;  % = 4*nnf, that is, 4 times the number of analysis frequencies.
df = 0.001                                  * time;

freqs = df*[[0:Nf/2-1] [-Nf/2:-1]].';

eval(['[s,a] = STASTurbine_'  nm ' ();']);
eval(['epar  = STASElectric_' nm ' ();']);
eval(['ppar  = STASPitch_'    nm ' ();']);
eval(['ypar  = STASYaw_'      nm ' ();']);
eval(['c     = STASControl_'  nm ' ();']);
eval(['m     = STASSensor_'   nm ' ();']);

a.dens = 1.225                     / density;
a.visc = 1.789e-5                  / viscosity;
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

blxdof = cell2mat (bldof(1));
bludof = cell2mat (bldof(2));
blydof = cell2mat (bldof(3));

Nx = size(xpsi,1);
Ny = size(ypsi,1);
Neta = 6 + s.foundation.Nmod + s.tower.Nmod + s.nacelle.Nmod ...
     + s.driveshaft.Nmod + s.blade(1).Nmod + s.blade(2).Nmod ...
     + s.blade(3).Nmod + 6;
Nret = size(dret,1);

iW = 2*Neta - 4;
WW = xpsi(iW);
a3P = 3*WW;

Nacp = size(a.icp,1);
Nxa = 7*3*Nacp;

% Get aero mode shapes.
qpsi   = ypsi(1:Ndj);
qdpsi  = ypsi(Ndj+[1:Ndj]);
qddpsi = ypsi(2*Ndj+[1:Ndj]);
Fpsi   = ypsi(3*Ndj+[1:Ndj]);
[b1,b2,b3] = MBCindices_Ndj (Ndj,idofs);
[TpsiB_Ndj,TBpsi_Ndj] = MBC (Ndj,b1,b2,b3,azi);
q   = TpsiB_Ndj*qpsi;
qd  = TpsiB_Ndj*qdpsi;
Psia = getPsia (0,s,a,azi,q,qd,P,shape0);

% Prepare the selected inputs: Pc from the plant control, rotor-
% average wind speed and direction, waterline wave force, and
% terminal voltage.
iuPc = Ndj + Neta + 3*Nae + 5;
inVx = Ndj + Neta + [1:3:3*Neb-2].';  % Collective components.
inVy = Ndj + Neta + [2:3:3*Neb-1].';
inwx = 6*(Nmud+Nwater-1) + 1;
inwy = 6*(Nmud+Nwater-1) + 2;
iuvs = Ndj + Neta + 3*Nae + [3:4].';

Nu = 6;

BB = sparse(Nx,Nu);
D0 = sparse(Ny,Nu);
BB(:,1) = Bpsi(:,iuPc);
D0(:,1) = Dpsi(:,iuPc);
mag = sqrt(Vg0(1)^2 + Vg0(2)^2);
ang = atan2(Vg0(2),Vg0(1));
ca = cos(ang);
sa = sin(ang);
BB(:,2:3) = [sum(Bpsi(:,inVx),2) sum(Bpsi(:,inVy),2)]*[ca, -mag*sa;sa, mag*ca];
D0(:,2:3) = [sum(Dpsi(:,inVx),2) sum(Dpsi(:,inVy),2)]*[ca, -mag*sa;sa, mag*ca];
mag = sqrt(Fw0(1)^2 + Fw0(2)^2);
ang = atan2(Fw0(2),Fw0(1));
ca = cos(ang);
sa = sin(ang);
BB(:,4) = [Bpsi(:,inwx) Bpsi(:,inwy)]*[ca;sa];
D0(:,4) = [Dpsi(:,inwx) Dpsi(:,inwy)]*[ca;sa];
BB(:,5:6) = Bpsi(:,iuvs);
D0(:,5:6) = Dpsi(:,iuvs);

% Prepare the selected outputs:
% Sensors: W, beta0, yaw, Pe, vnac, V, thV.
% Other outputs: a, FT, Pa, qfnd, qbl, ist. (Ignore Qa for the present study.)
ielf = Nmud + 1;
ifnd = idofs(1) + 6*(ielf-1) + [1:12].';

ielb = 1;
ibld = zeros(12,3);
ibld(:,1) = idofs(6) + 6*(ielb-1) + [1:12].';
ibld(:,2) = idofs(7) + 6*(ielb-1) + [1:12].';
ibld(:,3) = idofs(8) + 6*(ielb-1) + [1:12].';

Nyo = 61;
CCy = sparse(Nyo,Nx);
DDy = sparse(Nyo,Nu);

CCy(1,iW) = 1;                  % Rotor speed.
CCy(2,2*Neta+Nxa+1) = 1;        % Collective pitch.
CCy(3,2*Neta+Nxa+7) = 1;        % Yaw.
CCy(4,2*Neta+Nxa+8+25+15) = 1;  % Electric power.
CCy(5,2*Neta+Nxa+8+25+19) = 1;  % vnac,x.
CCy(6,2*Neta+Nxa+8+25+22) = 1;  % vnac,y.
DDy(7,2) = 1;                   % V.
DDy(8,3) = 1;                   % thV.

% da = (-1/V0) dVi + (Vi0/V0^2) dV*.
% This is complicated by the fact that the aero states are transformed.
% Transform the collective part back to element-by-element, and average
% over the swept area.
rbel = s.driveshaft.Lh + 0.5*a.Lel(1) ...
     + [0;cumsum(0.5*(a.Lel(1:Neb-1)+a.Lel(2:Neb)))];
r2 = rbel.^2;  % Weight according to swept area, proportional to r^2.
sr2 = sum(r2);
wgt = r2/sr2;
PsiVi = Psia(6:7:7*Neb-1,6:7:(Nxa/3)-1);
ind = 2*Neta + [6:7:(Nxa/3)-1].';
CCy(9,ind) = -sum(full(PsiVi).*wgt)/Vg0(1);
Vi0 = PsiVi*xpsi(ind);
Viavg = sum(Vi0.*wgt);
DDy(9,2) = Viavg/(Vg0(1)^2);

% Rotor thrust.
ind = 3*Ndj+[1:Ndj];
FB = TpsiB_Ndj*ypsi(ind);
CFB = TpsiB_Ndj*Cpsi(ind,:);    % Note, rows are now body, columns (states)
CqB = TpsiB_Ndj*Cpsi(1:Ndj,:);  % are still MBC.
DFB = TpsiB_Ndj*D0(ind,:);

FT = zeros(3*s.blade(1).Nnod,1);
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);
Td_n = [cp0 -sp0 0;sp0 cp0 0;0 0 1];
Try = Tn_y;
[Tyy0,dTyy0] = dTdth (q(idofs(3)+[4:6]));
Ty0g = TFromTheta (P(idofs(3)+[4:6]));
for ibl = 1:3

   jcn = s.blade(1).Nnod*(ibl-1);

   bldof = idofs(5+ibl);
   [Tpp0,dTpp0] = dTdth (q(bldof+[4:6]));
   Tp0g = TFromTheta (P(bldof+[4:6]));

   indF = bldof + [1:6*s.blade(1).Nnod].';
   indy = 3*Ndj + indF;

   TT = ((Ty0g*Tyy0*Try).')*Tp0g*Tpp0;
   for inod = 2:s.blade(1).Nnod
      ic6 = 6*(inod-1);
      ic3 = 3*(inod-1);
      Fr = TT*FB(indF(ic6+[1:3]));
      FT(jcn+inod) = Fr(3);

      CCy(10,:) = CCy(10,:)                                 ...
                + [0 0 1]*TT*CFB(indF(ic6+[1:3]),:);
      DDy(10,:) = DDy(10,:)                                 ...
                + [0 0 1]*TT*DFB(indF(ic6+[1:3]),:);

      for jj = 1:3
         jc3 = 3*(jj-1);
         CCy(10,:) = CCy(10,:)                                           ...
                   + [0 0 1]*((Ty0g*dTyy0(:,jc3+[1:3])*Try).')           ...
                   * Tp0g*Tpp0*FB(indF(ic6+[1:3]))*CqB(idofs(3)+3+jj,:)  ...
                   + [0 0 1]*((Ty0g*Tyy0*Try).')*Tp0g*dTpp0(:,jc3+[1:3]) ...
                   * FB(indF(ic6+[1:3]))*CqB(bldof+3+jj,:);
         % DqB = 0, so no contribution to DDy(10,:).
      end

   end

end

% Available power.
DDy(11,2) = CPstar*1.5*dens*Area*Vmag^2;

% Foundation and blade root q's.
CCy(12:23,:) = Cpsi(ifnd,:);
CCy(24:35,:) = Cpsi(ibld(:,1),:);
CCy(36:47,:) = Cpsi(ibld(:,2),:);
CCy(48:59,:) = Cpsi(ibld(:,3),:);

% Transformer terminal currents.
CCy(60:61,:) = Cpsi(4*Ndj+4+1+[1:2],:);

% Eliminate the azimuth DOF.
Nret = size(dret,1);

LA = Lpsi(dret,dret)\Apsi(dret,dret);
LB = Lpsi(dret,dret)\BB(dret,:);
CC = CCy(:,dret);
DD = DDy;

fkey = [0.01:0.01:fco].';
tol = 1.e-3;
Yflg = 1;
wrtflg = 1;
[Amod,Bmod,Cmod] = modalTransformation (LA,LB,CC,fkey,tol,Yflg,wrtflg);
Nxm = size(Amod,1);
Dmod = DD;

tol = 1.e-3;
fmatch = [0.01:0.01:fco].';
[Amr,Bmr,Cmr,Dmr] = modelReduction (Amod,Bmod,Cmod,Dmod,fmatch,tol);
%Amr = Amod;
%Bmr = Bmod;
%Cmr = Cmod;
%Dmr = Dmod;
Nxred = size(Amr,1);

eval(["save ('-binary','Amod" txt Vstr ".bin','Amod');"]);
eval(["save ('-binary','Bmod" txt Vstr ".bin','Bmod');"]);
eval(["save ('-binary','Cmod" txt Vstr ".bin','Cmod');"]);
eval(["save ('-binary','Dmod" txt Vstr ".bin','Dmod');"]);
eval(["save ('-binary','Amr" txt Vstr ".bin','Amr');"]);
eval(["save ('-binary','Bmr" txt Vstr ".bin','Bmr');"]);
eval(["save ('-binary','Cmr" txt Vstr ".bin','Cmr');"]);
eval(["save ('-binary','Dmr" txt Vstr ".bin','Dmr');"]);

% Check some TFs.
TFs = zeros(Nf/2,Nyo*Nu);
TFr = zeros(Nf/2,Nyo*Nu);
for ifreq = 1:Nf/2

   f = freqs(ifreq);
   w = 2*pi*f;

   TF = CC*((i*w*speye(Nret) - LA)\LB) + DD;
   TFs(ifreq,:) = reshape(TF,1,Nyo*Nu);

   TF = Cmr*((i*w*speye(Nxred) - Amr)\Bmr) + Dmr;
   TFr(ifreq,:) = reshape(TF,1,Nyo*Nu);

end

iu = 2;  % 1: Pc, 2: V, 3: thV, 4: Fw, 5-6 vs.
iy = 10;  % 1: W, 2: beta0, 3: yaw, 4: Pe, 5-6: vnac, 7: V, 8: thV,
         % 9: a, 10: FT, 11: Pa, 12-23: qfnd, 24-59: qbl.
icol = Nyo*(iu-1) + iy;

figure(1);
clf;
hold on;
semilogy(freqs(2:Nf/2),abs(TFs(2:Nf/2,icol)));
semilogy(freqs(2:Nf/2),abs(TFr(2:Nf/2,icol)));
hold off;

% Check eigenvalues of the reduced matrices.
[slap,shp,ifrq] = eigVal (Amr);

% Diagonalize and Y transform.
[AY,BY,CY,Phi,slap] = Ytransform (Amr,Bmr,Cmr);
DY = Dmr;

% Make high-frequency modes quasi-static.

ind = [1:Nxred].';
jqs = (abs(slap)/(2*pi) > fco);
iqs = ind(jqs);
[Ap,ret,jnk] = partitionMatrix (AY,iqs,iqs);
AA = AY(iqs,iqs)\AY(iqs,ret);
AB = AY(iqs,iqs)\BY(iqs,:);
Aq = AY(ret,ret) - AY(ret,iqs)*AA;
Bq = BY(ret,:) - AY(ret,iqs)*AB;
Cq = CY(:,ret) - CY(:,iqs)*AA;
Dq = DY - CY(:,iqs)*AB;

%Aq = AY;
%Bq = BY;
%Cq = CY;
%Dq = DY;
Nxq = size(Aq,1);

% Check some more TFs.
TFq = zeros(Nf/2,Nyo*Nu);
for ifreq = 1:Nf/2

   f = freqs(ifreq);
   w = 2*pi*f;

   TF = Cq*((i*w*speye(Nxq) - Aq)\Bq) + Dq;
   TFq(ifreq,:) = reshape(TF,1,Nyo*Nu);

end

figure(2);
clf;
hold on;
semilogy(freqs(2:Nf/2),abs(TFs(2:Nf/2,icol)));
semilogy(freqs(2:Nf/2),abs(TFr(2:Nf/2,icol)));
semilogy(freqs(2:Nf/2),abs(TFq(2:Nf/2,icol)));
hold off;

[AT,BT,dt] = SSDiscreteTime (2,Aq,Bq,dt);
[slap,shp,ifrq] = eigVal_silent (AT);
abs(slap)

eval(["save ('-binary','Aq" txt Vstr ".bin','Aq');"]);
eval(["save ('-binary','Bq" txt Vstr ".bin','Bq');"]);
eval(["save ('-binary','Cq" txt Vstr ".bin','Cq');"]);
eval(["save ('-binary','Dq" txt Vstr ".bin','Dq');"]);
eval(["save ('-binary','AT" txt Vstr ".bin','AT');"]);
eval(["save ('-binary','BT" txt Vstr ".bin','BT');"]);

%{

% Check some step functions.
Nt = 600;
yt = zeros(Nt,Nyo*Nu);
xt = zeros(Nxq,Nu);
u  = zeros(Nu,1);
ts = zeros(Nt,1);
for it = 2:Nt

   ts(it) = ts(it-1) + dt;

   if (it > 10)
      u = ones(Nu,1);
   end
   for iu = 1:Nu
      icn = Nyo*(iu-1);
      u1 = zeros(Nu,1);
      u1(iu) = u(iu);
      xt(:,iu) = AT*xt(:,iu) + BT*u1;
      yt(it,icn+[1:Nyo]) = (Cq*xt(:,iu) + Dq*u1).';
   end

end

% u(inputs) : 1: Pc, 2: V, 3: thV, 4: Fw, 5-6: vs
% y(sensors): 1: W, 2: beta0, 3: yaw, 4: Pe, 5-6: vnac, 7: V, 8: thV
% y(other)  : 9: a, 10: FT, 11: Pa, 12-23: qfnd, 24-59: qbl.
for iu = 1:4
   icn = Nyo*(iu-1);
   figure(iu+2);
   clf;
   for iyout = 1:8
      subplot(2,4,iyout);
      hold on;
      plot (ts,yt(:,icn+iyout),'color',[0.0, 0.0, 0.0],'linewidth',2.0);
      hold off;
   end
end

%}
