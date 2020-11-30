% Get spectra of internal loads at particular locations.

clear;

pkg load statistics;

[length,time,mass,current,voltage,        ...
 velocity,force,power,stress,ndens,nvisc, ...
 stiffness,damping,resistance,inductance, ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

load 'Sij_V10_I183_Lu180_Y0.bin';
load 'Sijp_V10_I183_Lu180_Y0.bin';
load 'cvg_V10_I183_Lu180_Y0.bin';
load 'Vavg_V10_I183_Lu180_Y0.bin';

load 'Savg_Hs20_TP6.bin';

nm = 'DTU10MW';
Vmag = 10;

yaw = 0;
betas = [0;0;0];
azi = 0;
thw = 0*pi/180;

Nmud = 5;

Nf = 2^12;  % = 4*nnf, that is, 4 times the number of analysis frequencies.
df = 0.001;

freqs = df*[[0:Nf/2-1] [-Nf/2:-1]].';

eval(['[s,a] = STASTurbine_'  nm ' ();']);
eval(['epar  = STASElectric_' nm ' ();']);
eval(['ppar  = STASPitch_'    nm ' ();']);
eval(['ypar  = STASYaw_'      nm ' ();']);
eval(['c     = STASControl_'  nm ' ();']);

[length,time,mass,current,voltage,        ...
 velocity,force,power,stress,ndens,nvisc, ...
 stiffness,damping,resistance,inductance, ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
[imdofs,Nmd] = getmdofRefs (s);
Ndj = Ndof + 6;

Pin = assemblePin (s);
[q,P,Ts0_B,TB0_g] =                                                     ...
      undeformedPosition (Pin,yaw,s.nacelle.delta,azi,s.driveshaft.phi, ...
                          betas,0,idofs,idofm,inods,inodm);

Nb  = a.Nb;
Neb = a.Neb;
Nae = Nb*Neb;

Nmud = s.foundation.Nmud;
Nwater = s.foundation.Nwater;

% Load the operating-point state vector and matrices.
txt = '_V';
eval(["load 'xpsi" txt int2str(10*Vmag) ".txt';"]);
eval(["load 'dxpsi" txt int2str(10*Vmag) ".txt';"]);
eval(["load 'ypsi" txt int2str(10*Vmag) ".txt';"]);
eval(["load 'Lpsi" txt int2str(10*Vmag) ".bin';"]);
eval(["load 'Apsi" txt int2str(10*Vmag) ".bin';"]);
eval(["load 'Bpsi" txt int2str(10*Vmag) ".bin';"]);
eval(["load 'Cpsi" txt int2str(10*Vmag) ".bin';"]);
eval(["load 'Dpsi" txt int2str(10*Vmag) ".bin';"]);
eval(["load 'dret" txt int2str(10*Vmag) ".bin';"]);
eval(["load 'bldof" txt int2str(round(10*Vmag)) ".bin';"]);
eval(["xpsi = xpsi" txt int2str(10*Vmag) ";"]);
eval(["dxpsi = dxpsi" txt int2str(10*Vmag) ";"]);
eval(["ypsi = ypsi" txt int2str(10*Vmag) ";"]);

blxdof = cell2mat (bldof(1));
bludof = cell2mat (bldof(2));
blydof = cell2mat (bldof(3));

Nx = size(xpsi,1);
Nret = size(dret,1);
Neta = 6 + s.foundation.Nmod + s.tower.Nmod + s.nacelle.Nmod ...
     + s.driveshaft.Nmod + s.blade(1).Nmod + s.blade(2).Nmod ...
     + s.blade(3).Nmod + 6;

WW = ypsi(Ndj+idofs(4)+6);

% Always partition/eliminate first, then invert L afterwards.
geteigs;

% Diagonalize.
Phi = shp(:,ifrq);
Psi = inv(Phi);
ii = [1:Nret].';
jj = ii;

% Some of the very closely spaced eigenfrequencies may lie out of order.
% This causes problems when I assume that conjugate modes are mirrored
% in the sorted mode shape matrix.  The solution is to force complex
% conjugacy by copying half the complex mode shapes and eigenvalues to
% the other half.  This is done in the 'for' loop below.
iY = speye(Nret);
Nosc = 0;
for im1 = 1:floor(Nret/2)

   imn = Nret - (im1-1);

   if (abs(imag(slap(im1))) > 0)
      Nosc = Nosc + 1;
      Phi(:,imn) = conj(Phi(:,im1));
      Psi(imn,:) = conj(Psi(im1,:));
      slap(imn) = conj(slap(im1));
      iY([im1 imn],[im1 imn]) = [1 1;-i i];
   end

end
Nexp = Nret - 2*Nosc;
Nmds = Nosc + Nexp;

Lam = sparse (ii,jj,slap,Nret,Nret);

% The wind is input as a set of rotationally-sampled velocity spectra
% in global coordinates, while the waves are input as unidirectional
% force spectra, in the wave direction, acting on the foundation nodes
% from the seabed to transition piece (z = surface+10 m for the DTU
% 10 MW).
%
% Prepare B matrices specifically for the wind and wave inputs.  The
% transform from wave coordinates to foundation body coordinates
% should be done as part of the Bw matrix.
Nwnod = s.foundation.Nnod - s.foundation.Nmud;
ct = cos(thw);
st = sin(thw);
Twg = [ct, -st,  0; ...
       st,  ct,  0; ...
        0,   0,  1];
Twgs = diagmat (2*Nwnod,Twg);
Twxg = Twgs(:,[1:6:6*Nwnod-5]);
LB = Lpsi(dret,dret)\Bpsi(dret,:);
Bv = LB(:,Ndj+Neta+[1:3*Nae]);
Bw = LB(:,idofs(1)+[6*s.foundation.Nmud+1:6*s.foundation.Nnod])*Twxg;

% I want to know the transfer functions between wind and wave inputs and
% internal loads; integrated wave loads and internal loads; control inputs
% and internal loads.  For internal loads, I have to convert from modal to
% nodal displacements, transform from body to section coordinates, and
% then use the element stiffness matrix to extract the internal loads.
[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

inac = idofs(3) + [1:2].';             % For nacelle motions.
iblt = idofs(6) + 6*Neb + [2:3].';     % For blade tip motions.

ielf = Nmud + 1;
ifnd = idofs(1) + 6*(ielf-1) + [1:12].';

ielt = 10;
itow = idofs(2) + 6*(ielt-1) + [1:12].';

ield = 4;
idrv = idofs(4) + 6*(ield-1) + [1:12].';

ielb = 1;
ibl = zeros(12,3);
ibl(:,1) = idofs(6) + 6*(ielb-1) + [1:12].';
ibl(:,2) = idofs(7) + 6*(ielb-1) + [1:12].';
ibl(:,3) = idofs(8) + 6*(ielb-1) + [1:12].';

% Foundation.
if (ielf == 1)
   qn1 = zeros(6,1);
   qn2 = ypsi(ifnd(7:12));
   Pn1 = [0;0;0;P(ifnd(10:12))];  % Same T_B^B0 orientation as second node.
   Pn2 = P(ifnd(7:12));
   kes = s.foundation.ke_s(:,12*(ielf-1)+[1:12]);
   conn = s.foundation.conn(:,ielf);
   [Pintf,dPintfdq] = internalLoads (qn1,qn2,Pn1,Pn2,kes);
   dPintfdq(:,1:6) = 0;
else
   qn1 = ypsi(ifnd(1:6));
   qn2 = ypsi(ifnd(7:12));
   Pn1 = P(ifnd(1:6));
   Pn2 = P(ifnd(7:12));
   kes = s.foundation.ke_s(:,12*(ielf-1)+[1:12]);
   conn = s.foundation.conn(:,ielf);
   [Pintf,dPintfdq] = internalLoads (qn1,qn2,Pn1,Pn2,kes);

end

% Tower.
if (ielt == 1)
   qn1 = zeros(6,1);
   qn2 = ypsi(itow(7:12));
   Pn1 = [0;0;0;P(itow(10:12))];  % Same T_B^B0 orientation as second node.
   Pn2 = P(itow(7:12));
   kes = s.tower.ke_s(:,12*(ielt-1)+[1:12]);
   conn = s.tower.conn(:,ielt);
   [Pintt,dPinttdq] = internalLoads (qn1,qn2,Pn1,Pn2,kes);
   dPinttdq(:,1:6) = 0;
else
   qn1 = ypsi(itow(1:6));
   qn2 = ypsi(itow(7:12));
   Pn1 = P(itow(1:6));
   Pn2 = P(itow(7:12));
   kes = s.tower.ke_s(:,12*(ielt-1)+[1:12]);
   conn = s.tower.conn(:,ielt);
   [Pintt,dPinttdq] = internalLoads (qn1,qn2,Pn1,Pn2,kes);
end

% Driveshaft.
if (ield == 1)
   qn1 = zeros(6,1);
   qn2 = ypsi(idrv(7:12));
   Pn1 = [0;0;0;P(idrv(10:12))];  % Same T_B^B0 orientation as second node.
   Pn2 = P(idrv(7:12));
   kes = s.driveshaft.ke_s(:,12*(ield-1)+[1:12]);
   conn = s.driveshaft.conn(:,ield);
   [Pintd,dPintddq] = internalLoads (qn1,qn2,Pn1,Pn2,kes);
   dPintddq(:,1:6) = 0;
else
   qn1 = ypsi(idrv(1:6));
   qn2 = ypsi(idrv(7:12));
   Pn1 = P(idrv(1:6));
   Pn2 = P(idrv(7:12));
   kes = s.driveshaft.ke_s(:,12*(ield-1)+[1:12]);
   conn = s.driveshaft.conn(:,ield);
   [Pintd,dPintddq] = internalLoads (qn1,qn2,Pn1,Pn2,kes);
end

% Blades.  These need to include the cross-spectra for the three blades.
% The qn's are based on ypsi, which is in MBC.  To appropriately match
% the P's I need to put things in body coordinates.
qpsi = ypsi(1:Ndj);
[b1,b2,b3] = MBCindices_Ndj (Ndj,idofs);
[TpsiB,TBpsi] = MBC (Ndj,b1,b2,b3,azi);
qB = TpsiB*qpsi;

Pintb = zeros(18,1);
dPintbdq = sparse(18,36);
qbs = zeros(36,1);
for ib = 1:3
   ir6 = 6*(ib-1);
   ic12 = 12*(ib-1);
   if (ielb == 1)
      qn1 = zeros(6,1);
      qn2 = qB(ibl(7:12,ib));
      Pn1 = [0;0;0;P(ibl(10:12,ib))];  % Same T_B^B0 orientation as second node.
      Pn2 = P(ibl(7:12,ib));
      kes = s.blade(ib).ke_s(:,12*(ielb-1)+[1:12]);
      conn = s.blade(ib).conn(:,ielb);
      [Pb,dPb] = internalLoads (qn1,qn2,Pn1,Pn2,kes);
      Pintb(ir6+[1:6]) = Pb;
      dPintbdq(ir6+[1:6],ic12+[1:12]) = dPb;
      dPintbdq(ir6+[1:6],ic12+[1:6]) = 0;
   else
      qn1 = qB(ibl(1:6,ib));
      qn2 = qB(ibl(7:12,ib));
      Pn1 = P(ibl(1:6,ib));
      Pn2 = P(ibl(7:12,ib));
      kes = s.blade(ib).ke_s(:,12*(ielb-1)+[1:12]);
      conn = s.blade(ib).conn(:,ielb);
      [Pintb,dPintbdq] = internalLoads (qn1,qn2,Pn1,Pn2,kes);
      Pintb(ir6+[1:6]) = Pb;
      dPintbdq(ir6+[1:6],ic12+[1:12]) = dPb;
   end
   qbs(ic12+[1:12]) = [qn1;qn2];  % Save for later.
end

% Transform the mean loads and disps to MBC, needed to get the right 1P.
b1 = [1:6].';
b2 = [7:12].';
b3 = [13:18].';
[TpsiB,TBpsi] = MBC (18,b1,b2,b3,azi);
Pbpsi = TBpsi*Pintb;
b1 = [1:12].';
b2 = [13:24].';
b3 = [25:36].';
[TpsiB,TBpsi] = MBC (36,b1,b2,b3,azi);
qbpsi = TBpsi*qbs;

PBv = Psi*Bv;
PBw = Psi*Bw;

% I want internal load spectra in the blades and driveshaft in rotating
% coordinates.  The dPintbdq transfer function is expressed in body 
% coordinates.  The dPintddq transfer function should be expressed in
% the coordinates associated with the current rotor position.
%
% So, apply the MBC transfer functions to get the required qpsi cross-
% spectra.  Then transform these Sqpsi spectra to rotating coordinates,
% for the driveshaft and blades.  Then apply the internal-load TFs.
qdofs = [ifnd;itow;idrv;ibl(:,1);ibl(:,2);ibl(:,3)];
CC = Cpsi(qdofs,dret);
CP = CC*Phi;

% Displacement spectra.
%
% 12-by-12 cross-spectra for each of three components (foundation,
% tower, driveshaft), and 36-by-36 for the blades.
Nyc = 12*3 + 36;
Nys = 144*3 + 36^2; 

i3a = [1:3:3*Nae-2].';
i3b = [2:3:3*Nae-1].';
i3c = [3:3:3*Nae].';

i9a = [1:9:9*(Nae^2)-8].';
i9b = [2:9:9*(Nae^2)-7].';
i9c = [3:9:9*(Nae^2)-6].';
i9d = [4:9:9*(Nae^2)-5].';
i9e = [5:9:9*(Nae^2)-4].';
i9f = [6:9:9*(Nae^2)-3].';
i9g = [7:9:9*(Nae^2)-2].';
i9h = [8:9:9*(Nae^2)-1].';
i9i = [9:9:9*(Nae^2)].';

Sq = zeros(Nf,Nys);
cq = zeros(Nf,Nyc);
fidm = fopen("TFmag.txt","w");
fidp = fopen("TFph.txt","w");
for ifreq = 1:Nf

%printf('%6d of %6d\n',ifreq,Nf);
%fflush(stdout);

   f = freqs(ifreq);
   w = 2*pi*f;

   iwslap = i*w*speye(Nret) - Lam;

   Hv = CP*(iwslap\PBv);
   Hw = CP*(iwslap\PBw);

valry = Hv(36+11,3*11+1);
valrz = Hv(36+12,3*11+1);
fprintf(fidm,"%+5.6e %+5.6e %+5.6e\n",f,abs(valry),abs(valrz));
fprintf(fidp,"%+5.6e %+5.6e %+5.6e\n",f, ...
        atan2(imag(valry),real(valry))/pi, ...
        atan2(imag(valrz),real(valrz))/pi);


%   Sijf = zeros(3*Nae,3*Nae);
%   for iel1 = 1:Nae
%      ir3 = 3*(iel1-1);
%      iref = 9*Nae*(iel1-1);
%      for iel2 = 1:Nae
%         ic3 = 3*(iel2-1);
%         icol = iref + 9*(iel2-1);
%         Sijf(ir3+[1:3],ic3+[1:3]) = reshape (Sij(ifreq,icol+[1:9]),3,3);
%      end
%   end

   % Gives the same answer as above, 10 times faster.
   Sijf = zeros(3*Nae,3*Nae);
   Sijf(i3a,i3a) = reshape(Sij(ifreq,i9a),Nae,Nae).';
   Sijf(i3b,i3a) = reshape(Sij(ifreq,i9b),Nae,Nae).';
   Sijf(i3c,i3a) = reshape(Sij(ifreq,i9c),Nae,Nae).';
   Sijf(i3a,i3b) = reshape(Sij(ifreq,i9d),Nae,Nae).';
   Sijf(i3b,i3b) = reshape(Sij(ifreq,i9e),Nae,Nae).';
   Sijf(i3c,i3b) = reshape(Sij(ifreq,i9f),Nae,Nae).';
   Sijf(i3a,i3c) = reshape(Sij(ifreq,i9g),Nae,Nae).';
   Sijf(i3b,i3c) = reshape(Sij(ifreq,i9h),Nae,Nae).';
   Sijf(i3c,i3c) = reshape(Sij(ifreq,i9i),Nae,Nae).';
  
 
   Swf = reshape (Savg(ifreq,:),Nwnod,Nwnod).'/(force^2/time);

   for icomp = 1:3

      ir12  = 12*(icomp-1);
      ir144 = 144*(icomp-1);

      Sqv = Hv(ir12+[1:12],:)*Sijf*(Hv(ir12+[1:12],:)');
      Sqw = Hw(ir12+[1:12],:)*Swf*(Hw(ir12+[1:12],:)');
      Sq(ifreq,ir144+[1:144]) = reshape(Sqv + Sqw,1,144);

      cq(ifreq,ir12+[1:12]) = (Hv(ir12+[1:12],:)*cvg(ifreq,:).').';

   end

   ir12  = 3*12;
   ir144 = 3*144;
   Sqv = Hv(ir12+[1:36],:)*Sijf*(Hv(ir12+[1:36],:)');
   Sqw = Hw(ir12+[1:36],:)*Swf*(Hw(ir12+[1:36],:)');
   Sq(ifreq,ir144+[1:36^2]) = reshape(Sqv + Sqw,1,36^2);   

   cq(ifreq,ir12+[1:36]) = (Hv(ir12+[1:36],:)*cvg(ifreq,:).').';

end
fclose('all');

% First convert the blade Sq's to rotating coordinates, retaining only
% Blade 1.
Sqbr = zeros(Nf,144);
for icol = 1:12
   for irow = 1:12
      icomp = 12*(icol-1) + irow;
      ind = [     irow + 36*(icol-1) + [0;36*12;36*24]; ...
             12 + irow + 36*(icol-1) + [0;36*12;36*24]; ...
             24 + irow + 36*(icol-1) + [0;36*12;36*24]];
      SS = Sq(:,3*144+ind);
      Stemp = MBCtoBodySpectra(1,SS,azi,WW,df);
      Sqbr(:,icomp) = Stemp(:,1);
   end
end




save('cq.bin','cq');
save('qbpsi.bin','qbpsi');




cqbr = sparse(Nf,12);
for icomp = 1:12
   ind = icomp + [0;12;24];
   [ctemp,qtemp] = MBCtoBodyCoeffs (1,cq(:,36+ind),qbpsi(ind),azi,WW,df);
   cqbr(:,icomp) = ctemp(:,1);
end

Spf = zeros(Nf,36);
Spt = zeros(Nf,36);
Spd = zeros(Nf,36);
Spb = zeros(Nf,36);
cpb = zeros(Nf,6);
for ifreq = 1:Nf

   Stemp = reshape(Sq(ifreq,1:144),12,12);
   Spf(ifreq,:) = reshape(dPintfdq*Stemp*(dPintfdq'),1,36);

   Stemp = reshape(Sq(ifreq,144+[1:144]),12,12);
   Spt(ifreq,:) = reshape(dPinttdq*Stemp*(dPinttdq'),1,36);

   Stemp = reshape(Sq(ifreq,2*144+[1:144]),12,12);
   Spd(ifreq,:) = reshape(dPintddq*Stemp*(dPintddq'),1,36);

   Stemp = reshape(Sqbr(ifreq,[1:144]),12,12);
   Spb(ifreq,:) = reshape(dPintbdq(1:6,1:12)*Stemp*(dPintbdq(1:6,1:12)'),1,36);

   cpb(ifreq,:) = (dPintbdq(1:6,1:12)*cqbr(ifreq,:).').';
   Stemp = reshape((cpb(ifreq,:).')*conj(cpb(ifreq,:))/df,1,36);
   Spb(ifreq,:) = Spb(ifreq,:) + Stemp;
%   Spb(ifreq,:) = Spb(ifreq,:) + abs(Stemp);

end

% Transform the driveshaft load spectra.  (Can't do this transformation
% on q because the displacement DOF's don't transform as vectors.  Loads,
% on the other hand, do.)
Spd0 = Spd;
for ii = 1:2
   ir3 = 3*(ii-1);
   for jj = 1:2
      ic3 = 3*(jj-1);
      ind = 6*ic3 + ir3 + [1 2 3 7 8 9 13 14 15].';
      SS = Spd0(:,ind);
      Spd(:,ind) = MBCtoBodySpectra(2,SS,azi,WW,df);

%'------------------------------'
%ifreq = 5;
%[ii jj]
%reshape(SS(5,:),3,3).'
%reshape(Spd(5,ind),3,3).'

   end
end

ind = [Nf/2+[1:Nf/2] [1:Nf/2]].';

iel1 = 10;
iel2 = 10;
iref = 9*Nae*(iel1-1) + 9*(iel2-1);
figure(1);
clf;
hold on;
semilogy(freqs(ind),abs(Sij(ind,iref+[1 5 9])));
hold off;

figure(2);
clf;
hold on;
semilogy(freqs(ind),abs(Savg(ind,9*16+10)));
hold off;

figure(3);
clf;
hold on;
semilogy(freqs(ind),abs(Sq(ind,1)));
semilogy(freqs(ind),abs(Sq(ind,13)));
semilogy(freqs(ind),abs(Sq(ind,3*144+2*36+3)));
hold off;

figure(4);
clf;
hold on;
semilogy(freqs(ind),abs(Spf(ind,6*3+4)));
semilogy(freqs(ind),abs(Spf(ind,6*4+5)));
semilogy(freqs(ind),abs(Spf(ind,6*5+6)));
hold off;

figure(5);
clf;
hold on;
semilogy(freqs(ind),abs(Spt(ind,6*3+4)));
semilogy(freqs(ind),abs(Spt(ind,6*4+5)));
semilogy(freqs(ind),abs(Spt(ind,6*5+6)));
hold off;

figure(6);
clf;
hold on;
semilogy(freqs(ind),abs(Spd(ind,6*3+4)));
semilogy(freqs(ind),abs(Spd(ind,6*4+5)));
semilogy(freqs(ind),abs(Spd(ind,6*5+6)));
hold off;

figure(7);
clf;
hold on;
semilogy(freqs(ind),abs(Spb(ind,6*3+4)));
semilogy(freqs(ind),abs(Spb(ind,6*4+5)));
semilogy(freqs(ind),abs(Spb(ind,6*5+6)));
hold off;

%{

fida = fopen('specsa.txt','w');
fidw = fopen('specsw.txt','w');
fidq = fopen('specsq.txt','w');
fidf = fopen('specsf.txt','w');
fidt = fopen('specst.txt','w');
fidd = fopen('specsd.txt','w');
for jj = 1:Nf

   if (jj <= Nf/2)
      ifreq = Nf/2 + jj;
   else
      ifreq = jj - Nf/2;
   end

   f = freqs(ifreq);

   fprintf(fid,'%+5.6e %+5.6e %+5.6e %+5.6e %+5.6e\n', ...
           Sij(),);

end

%}


% Stresses.

