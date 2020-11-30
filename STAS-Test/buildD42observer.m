
clear;

pkg load statistics;
pkg load control;

[length,time,mass,current,voltage,              ...
 velocity,force,power,stress,density,viscosity, ...
 stiffness,damping,resistance,inductance,       ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

%load 'Sij_V10_I183_Lu180_Y0.bin';
%load 'Sijp_V10_I183_Lu180_Y0.bin';
%load 'cvg_V10_I183_Lu180_Y0.bin';
%load 'Vavg_V10_I183_Lu180_Y0.bin';

%load 'Savg_Hs20_TP6.bin';

nm = 'DTU10MW';
inpnm = '_P060_';

modalAero = 0;

Vmag = 10                                   / (length/time);
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
atw = 0.01*(2*pi);
Fw0 = 1e6*[cos(wang), sin(wang)].'          / force;

aV = 0.05*(2*pi)                            * time;
% a3P computed based on W.
zet3P = 0.1;

fco = 0.5                                   * time;  % > fco solved as static.
dt = 0.05                                   / time;  % Set to -1 for auto.

Nsp = 7; % Get stresses at Nsp points over a quadrant of the circumference.
ths = (pi/12)*[0:Nsp-1].';
ct = cos(ths);
st = sin(ths);

Nf = 2^12;  % = 4*nnf, that is, 4 times the number of analysis frequencies.
df = 0.001                                  * time;

freqs = df*[[0:Nf/2-1] [-Nf/2:-1]].';

eval(['[s,a] = STASTurbine_'  nm ' ();']);
eval(['epar  = STASElectric_' nm ' ();']);
eval(['ppar  = STASPitch_'    nm ' ();']);
eval(['ypar  = STASYaw_'      nm ' ();']);
eval(['c     = STASControl_'  nm ' ();']);

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
[imdofs,Nmd] = getmdofRefs (s);
Ndj = Ndof + 6;
Nnod = Ndof/6;

Pin = assemblePin (s);
[jnkq,P,Ts0_B,TB0_g] =                                                  ...
      undeformedPosition (Pin,yaw,s.nacelle.delta,azi,s.driveshaft.phi, ...
                          betas,0,idofs,idofm,inods,inodm);

Nb  = a.Nb;
Neb = a.Neb;
Nae = Nb*Neb;

Nmud = s.foundation.Nmud;
Nwater = s.foundation.Nwater;

% Load the operating-point state vector and matrices.
txt = inpnm;
if (Vmag >= 10)
   Vstr = ['V' int2str(round(10*Vmag))];
else
   Vstr = ['V0' int2str(round(10*Vmag))];
end
eval(["load 'xop" txt Vstr ".txt';"]);
eval(["load 'yop" txt Vstr ".txt';"]);
eval(["load 'Lop" txt Vstr ".bin';"]);
eval(["load 'Aop" txt Vstr ".bin';"]);
eval(["load 'Bop" txt Vstr ".bin';"]);
eval(["load 'Cop" txt Vstr ".bin';"]);
eval(["load 'Dop" txt Vstr ".bin';"]);
eval(["load 'shape" txt Vstr ".bin';"]);
eval(["load 'bldof" txt Vstr ".bin';"]);
eval(["xop = xop" txt Vstr ";"]);
eval(["yop = yop" txt Vstr ";"]);

blxdof = cell2mat (bldof(1));
bludof = cell2mat (bldof(2));
blydof = cell2mat (bldof(3));

Nx = size(xop,1);
Ny = size(yop,1);
Neta = 6 + s.foundation.Nmod + s.tower.Nmod + s.nacelle.Nmod ...
     + s.driveshaft.Nmod + s.blade(1).Nmod + s.blade(2).Nmod ...
     + s.blade(3).Nmod + 6;

iW = 2*Neta - 4;
WW = xop(iW);
a3P = 3*WW;

if (modalAero == 1)
   a.icp = [-1;-3];
end
Nacp = size(a.icp,1);
Nxa = 7*3*Nacp;

% Get aero mode shapes.
qpsi   = yop(1:Ndj);
qdpsi  = yop(Ndj+[1:Ndj]);
qddpsi = yop(2*Ndj+[1:Ndj]);
Fpsi   = yop(3*Ndj+[1:Ndj]);
[b1,b2,b3] = MBCindices_Ndj (Ndj,idofs);
[TpsiB_Ndj,TBpsi_Ndj] = MBC (Ndj,b1,b2,b3,azi);
q   = TpsiB_Ndj*qpsi;
qd  = TpsiB_Ndj*qdpsi;
Psia = getPsia (modalAero,s,a,azi,q,qd,P,shape0);

% Prepare the observer inputs: P, beta, yaw from the plant control, 
% rotor-average wind speed and direction, and waterline wave force.
iuP  = Ndj + Neta + 3*Nae + 4 + 8 + 1;
iub  = Ndj + Neta + 3*Nae + 1;
iuy  = Ndj + Neta + 3*Nae + 4;
inVx = Ndj + Neta + [1:3:3*Neb-2].';  % Collective components.
inVy = Ndj + Neta + [2:3:3*Neb-1].';
inwx = 6*(Nmud+Nwater-1) + 1;
inwy = 6*(Nmud+Nwater-1) + 2;

Nu = 6;

BB = sparse(Nx,Nu);
DDy = sparse(Ny,Nu);
BB (:,1) = Bop(:,iuP);
DDy(:,1) = Dop(:,iuP);
BB (:,2) = Bop(:,iub);
DDy(:,2) = Dop(:,iub);
BB (:,3) = Bop(:,iuy);
DDy(:,3) = Dop(:,iuy);
mag = sqrt(Vg0(1)^2 + Vg0(2)^2);
ang = atan2(Vg0(2),Vg0(1));
ca = cos(ang);
sa = sin(ang);
BB(:,4:5) = [sum(Bop(:,inVx),2) sum(Bop(:,inVy),2)]*[ca, -mag*sa;sa, mag*ca];
DDy(:,4:5) = [sum(Dop(:,inVx),2) sum(Dop(:,inVy),2)]*[ca, -mag*sa;sa, mag*ca];
mag = sqrt(Fw0(1)^2 + Fw0(2)^2);
ang = atan2(Fw0(2),Fw0(1));
ca = cos(ang);
sa = sin(ang);
BB(:,6) = [Bop(:,inwx) Bop(:,inwy)]*[ca;sa];
DDy(:,6) = [Dop(:,inwx) Dop(:,inwy)]*[ca;sa];

% Augment the matrices with the environmental models.
a32 = a3P^2;
tza = 2*zet3P*a3P;
aw2 = aw^2;
tzaw = 2*zetw*aw;
aae = [-aV, 0,    0,    0,    0,    0; ...
       0, -aV,    0,    0,    0,    0; ...
       0,   0,    0,    1,    0,    0; ...
       0,   0, -a32, -tza,    0,    0; ...
       0,   0,    0,    0,    0,    1; ...
       0,   0,    0,    0, -aw2, -tzaw];
bbe = [aV,  0,  0; ...
       0,  aV,  0; ...
       0,   0,  0; ...
       tza, 0,  0; ...
       0,   0,  0; ...
       0,   0,  tzaw];

Nxe = size(aae,1);
Nxaug = Nx + Nxe;
Laug = sparse(Nxaug,Nxaug);
Aaug = sparse(Nxaug,Nxaug);
Baug = sparse(Nxaug,Nu);
Caug = sparse(Ny,Nxaug);
Aaug(1:Nx,1:Nx) = Aop;
Laug(1:Nx,1:Nx) = Lop;
Baug(1:Nx,1:3) = BB(:,1:3);         % Only the control columns.
Caug(:,1:Nx) = Cop;
Aaug(Nx+[1:Nxe],Nx+[1:Nxe]) = aae;
Baug(Nx+[1:Nxe],4:6) = bbe;         % External V,thV,Fw feed the env model.
Aaug(1:Nx,Nx+1) = BB(:,4);          % Env model x feeds the rest of the model.
Aaug(1:Nx,Nx+2) = BB(:,5);
%Aaug(1:Nx,Nx+4) = BB(:,4);
Aaug(1:Nx,Nx+6) = BB(:,6);
Caug(:,Nx+1) = DDy(:,4);
Caug(:,Nx+2) = DDy(:,5);
%Caug(:,Nx+4) = DDy(:,4);
Caug(:,Nx+6) = DDy(:,6);
Laug(Nx+[1:Nxe],Nx+[1:Nxe]) = speye(Nxe);

% Prepare the observer outputs:
% Sensors: W, beta0, yaw, Pe, vnac, V, thV.
% Other outputs: a, FT, Pa, qfnd, qbl. (Ignore Qa for the present study.)
ielf = Nmud + 1;
ifnd = idofs(1) + 6*(ielf-1) + [1:12].';

ielb = 1;
ibld = zeros(12,3);
ibld(:,1) = idofs(6) + 6*(ielb-1) + [1:12].';
ibld(:,2) = idofs(7) + 6*(ielb-1) + [1:12].';
ibld(:,3) = idofs(8) + 6*(ielb-1) + [1:12].';

Nyo = 61;
ccy = sparse(Nyo,Nxaug);

ind = 2*Neta + Nxa + 8 + 25 + [1, 2, 5, 6, 7, 8, 9, 10].';
ccy(1:8,ind) = speye(8);

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
ccy(9,ind) = -sum(full(PsiVi).*wgt)/Vg0(1);
Vi0 = PsiVi*xop(ind);
Viavg = sum(Vi0.*wgt);
ind = Nx + 1; % V*
ccy(9,ind) = Viavg/(Vg0(1)^2);

% Rotor thrust.
ind = 3*Ndj+[1:Ndj];
FB = TpsiB_Ndj*yop(ind);
CFB = TpsiB_Ndj*Caug(ind,:);    % Note, rows are now body, columns (states)
CqB = TpsiB_Ndj*Caug(1:Ndj,:);  % are still MBC.
DFB = TpsiB_Ndj*DDy(ind,:);

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

      ccy(10,:) = ccy(10,:)                                 ...
                + [0 0 1]*TT*CFB(indF(ic6+[1:3]),:);
% Is this double-counting?  Already put DDy's in Caug.
%      ccy(10,Nx+1) = ccy(10,Nx+1)                           ...  % V
%                   + [0 0 1]*TT*DFB(indF(ic6+[1:3]),4);
%      ccy(10,Nx+4) = ccy(10,Nx+4)                           ...  % V3P
%                + [0 0 1]*TT*DFB(indF(ic6+[1:3]),4);

      for jj = 1:3
         jc3 = 3*(jj-1);
         ccy(10,:) = ccy(10,:)                                           ...
                   + [0 0 1]*((Ty0g*dTyy0(:,jc3+[1:3])*Try).')           ...
                   * Tp0g*Tpp0*FB(indF(ic6+[1:3]))*CqB(idofs(3)+3+jj,:)  ...
                   + [0 0 1]*((Ty0g*Tyy0*Try).')*Tp0g*dTpp0(:,jc3+[1:3]) ...
                   * FB(indF(ic6+[1:3]))*CqB(bldof+3+jj,:);
         % DqB = 0, so no additional contribution.
      end

   end

end

% Available power.
ind = 2*Neta + Nxa + 8 + 25 + 9;
ccy(11,ind) = CPstar*1.5*dens*Area*(xop(ind)^2);

% Foundation and blade root q's.
ccy(12:23,:) = Caug(ifnd,:);
ccy(24:35,:) = Caug(ibld(:,1),:);
ccy(36:47,:) = Caug(ibld(:,2),:);
ccy(48:59,:) = Caug(ibld(:,3),:);

% Transformer terminal currents.
ccy(60:61,:) = Caug(4*Ndj+4+1+[1:2],:);

% Eliminate the azimuth DOF.
dret = [[1:Neta-6] Neta-[5 3 2 1 0] [(Neta+1):Nxaug]].';
Nret = size(dret,1);

LA = Laug(dret,dret)\Aaug(dret,dret);
LB = Laug(dret,dret)\Baug(dret,:);
CC = ccy(:,dret);

fkey = [0.01:0.01:fco].';
tol = 1.e-4;
Yflg = 1;
wrtflg = 1;
[Aobs,Bobs,Cobs] = modalTransformation (LA,LB,CC,fkey,tol,Yflg,wrtflg);
%Aobs = LA;
%Bobs = LB;
%Cobs = CC;
Dobs = sparse(size(Cobs,1),size(Bobs,2));
Nxm = size(Aobs,1);

tol = 1.e-3;
fmatch = [0.01:0.01:fco].';
[Ared,Bred,Cred,Dred] = modelReduction (Aobs,Bobs,Cobs,Dobs,fmatch,tol);
%Ared = Aobs;
%Bred = Bobs;
%Cred = Cobs;
%Dred = Dobs;
Nxred = size(Ared,1);

eval(["save ('-binary','Aobs" txt Vstr ".bin','Aobs');"]);
eval(["save ('-binary','Bobs" txt Vstr ".bin','Bobs');"]);
eval(["save ('-binary','Cobs" txt Vstr ".bin','Cobs');"]);
eval(["save ('-binary','Dobs" txt Vstr ".bin','Dobs');"]);
eval(["save ('-binary','Ared" txt Vstr ".bin','Ared');"]);
eval(["save ('-binary','Bred" txt Vstr ".bin','Bred');"]);
eval(["save ('-binary','Cred" txt Vstr ".bin','Cred');"]);
eval(["save ('-binary','Dred" txt Vstr ".bin','Dred');"]);

%{

% Check some TFs.
TFs = zeros(Nf/2,Nyo*Nu);
TFr = zeros(Nf/2,Nyo*Nu);
for ifreq = 1:Nf/2

   f = freqs(ifreq);
   w = 2*pi*f;

%   TF = Cobs*((i*w*speye(Nxm) - Aobs)\Bobs) + Dobs;
%   TFs(ifreq,:) = reshape(TF,1,Nyo*Nu);

   TF = CC*((i*w*speye(Nret) - LA)\LB);
   TFs(ifreq,:) = reshape(TF,1,Nyo*Nu);

   TF = Cred*((i*w*speye(Nxred) - Ared)\Bred) + Dred;
   TFr(ifreq,:) = reshape(TF,1,Nyo*Nu);

end

% 1: Pc, 2: beta, 3: yaw, 4: V, 5: thV, 6: Fw.
% 1: W, 2: beta0, 3: yaw, 4: Pe, 5-6: vnac, 7: V, 8: thV,
% 9: a, 10: FT, 11: Pa, 12-23: qfnd, 24-59: qbl.
for iu = 1:6

   icn = Nyo*(iu-1);

   figure(iu);
   clf;
   subplot(2,2,1)
   hold on;
   semilogy(freqs(2:Nf/2),abs(TFs(2:Nf/2,icn+[1, 2, 3])), ...
            'color',[0, 0, 0],'linewidth',3.0);
   semilogy(freqs(2:Nf/2),abs(TFr(2:Nf/2,icn+[1, 2, 3])), ...
            'color',[0.6, 0.6, 0.6],'linewidth',2.0);
   hold off;
   subplot(2,2,2)
   hold on;
   semilogy(freqs(2:Nf/2),abs(TFs(2:Nf/2,icn+[4, 9, 10, 11])), ...
            'color',[0, 0, 0],'linewidth',3.0);
   semilogy(freqs(2:Nf/2),abs(TFr(2:Nf/2,icn+[4, 9, 10, 11])), ...
            'color',[0.6, 0.6, 0.6],'linewidth',2.0);
   hold off;
   subplot(2,2,3)
   hold on;
   semilogy(freqs(2:Nf/2),abs(TFs(2:Nf/2,icn+[5:8])), ...
            'color',[0, 0, 0],'linewidth',3.0);
   semilogy(freqs(2:Nf/2),abs(TFr(2:Nf/2,icn+[5:8])), ...
            'color',[0.6, 0.6, 0.6],'linewidth',2.0);
   hold off;
   subplot(2,2,4)
   hold on;
   semilogy(freqs(2:Nf/2),abs(TFs(2:Nf/2,icn+[36:41])), ...
            'color',[0, 0, 0],'linewidth',3.0);
   semilogy(freqs(2:Nf/2),abs(TFr(2:Nf/2,icn+[36:41])), ...
            'color',[0.6, 0.6, 0.6],'linewidth',2.0);
   hold off;

end

%}

% Check eigenvalues of the reduced matrices.
[slap,shp,ifrq] = eigVal (Ared);

% Diagonalize and Y transform.
[AY,BY,CY,Phi,slap] = Ytransform (Ared,Bred,Cred);
DY = Dred;

% Make high-frequency modes quasi-static.

ind = [1:Nxred].';
jqs = (abs(slap)/(2*pi) > fco);
iqs = ind(jqs);
[Ap,ret,jnk] = partitionMatrix (AY,iqs,iqs);
AA = AY(iqs,iqs)\AY(iqs,ret);
AB = AY(iqs,iqs)\BY(iqs,:);
Ao = AY(ret,ret) - AY(ret,iqs)*AA;
Bo = BY(ret,:) - AY(ret,iqs)*AB;
Co = CY(:,ret) - CY(:,iqs)*AA;
Do = DY - CY(:,iqs)*AB;

%Ao = AY;
%Bo = BY;
%Co = CY;
%Do = DY;
Nxo = size(Ao,1);

%{

% Check some more TFs.
TFq = zeros(Nf/2,Nyo*Nu);
for ifreq = 1:Nf/2

   f = freqs(ifreq);
   w = 2*pi*f;

   TF = Co*((i*w*speye(Nxo) - Ao)\Bo) + Do;
   TFq(ifreq,:) = reshape(TF,1,Nyo*Nu);

end

for iu = 1:6

   icn = Nyo*(iu-1);

   figure(iu+6);
   clf;
   subplot(2,2,1)
   hold on;
   semilogy(freqs(2:Nf/2),abs(TFs(2:Nf/2,icn+[1, 2, 3])), ...
            'color',[0, 0, 0],'linewidth',3.0);
   semilogy(freqs(2:Nf/2),abs(TFq(2:Nf/2,icn+[1, 2, 3])), ...
            'color',[0.6, 0.6, 0.6],'linewidth',2.0);
   hold off;
   subplot(2,2,2)
   hold on;
   semilogy(freqs(2:Nf/2),abs(TFs(2:Nf/2,icn+[4, 9, 10, 11])), ...
            'color',[0, 0, 0],'linewidth',3.0);
   semilogy(freqs(2:Nf/2),abs(TFq(2:Nf/2,icn+[4, 9, 10, 11])), ...
            'color',[0.6, 0.6, 0.6],'linewidth',2.0);
   hold off;
   subplot(2,2,3)
   hold on;
   semilogy(freqs(2:Nf/2),abs(TFs(2:Nf/2,icn+[5:8])), ...
            'color',[0, 0, 0],'linewidth',3.0);
   semilogy(freqs(2:Nf/2),abs(TFq(2:Nf/2,icn+[5:8])), ...
            'color',[0.6, 0.6, 0.6],'linewidth',2.0);
   hold off;
   subplot(2,2,4)
   hold on;
   semilogy(freqs(2:Nf/2),abs(TFs(2:Nf/2,icn+[36:41])), ...
            'color',[0, 0, 0],'linewidth',3.0);
   semilogy(freqs(2:Nf/2),abs(TFq(2:Nf/2,icn+[36:41])), ...
            'color',[0.6, 0.6, 0.6],'linewidth',2.0);
   hold off;

end

%}

[ADT,BDT,dt] = SSDiscreteTime (2,Ao,Bo,dt);
[slap,shp,ifrq] = eigVal_silent (ADT);
abs(slap)

eval(["save ('-binary','Ao" txt Vstr ".bin','Ao');"]);
eval(["save ('-binary','Bo" txt Vstr ".bin','Bo');"]);
eval(["save ('-binary','Co" txt Vstr ".bin','Co');"]);
eval(["save ('-binary','Do" txt Vstr ".bin','Do');"]);
eval(["save ('-binary','ADT" txt Vstr ".bin','ADT');"]);
eval(["save ('-binary','BDT" txt Vstr ".bin','BDT');"]);

%{

% Check some step functions.
Nt = 600;
yo = zeros(Nt,Nyo*Nu);
xo = zeros(Nxo,Nu);
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
      xo(:,iu) = ADT*xo(:,iu) + BDT*u1;
      yo(it,icn+[1:Nyo]) = (Co*xo(:,iu) + Do*u1).';
   end

end

% 1: Pc, 2: beta, 3: yaw, 4: V, 5: thV, 6: Fw.
% y(sensors): 1: W, 2: beta0, 3: yaw, 4: Pe, 5-6: vnac, 7: V, 8: thV
% y(other)  : 9: a, 10: FT, 11: Pa, 12-23: qfnd, 24-59: qbl.
for iu = 1:6
   icn = Nyo*(iu-1);
   figure(iu+12);
   clf;
   for iyout = 1:8
      subplot(2,4,iyout);
      hold on;
      plot (ts,yo(:,icn+iyout),'color',[0.0, 0.0, 0.0],'linewidth',2.0);
      hold off;
   end
end

%}