% Generate turbulence spectra, including the effects of periodic
% wind shear and tower shadow.
%
% Missing functionality: (See previous STAS Wind versions.)
% - Tower-blade cross-spectra.
% - Real/Imaginary spectra for complex-step gradients.
% - The spectral scaling method for Vinf, W.

clear;

nm = 'DTU10MW';
eval(['[s,a] = STASTurbine_'  nm ' ();']);
load 'xpsi_P060_V100.txt';
load 'ypsi_P060_V100.txt';
xpsi = xpsi_P060_V100;
ypsi = ypsi_P060_V100;

[length,time,mass,current,voltage,        ...
 velocity,force,power,stress,ndens,nvisc, ...
 stiffness,damping,resistance,inductance, ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

%-------------------------------------------------------------------
% Inputs:

MBCflag = 1;

Vmag = 10.;
yaw = 0*pi/180;
Vinf = [Vmag*cos(yaw);Vmag*sin(yaw);0];
betas = [0;0;0];
psi0 = 0;

%TI    = 0.14*(0.75*Vmag + 5.6)./Vmag; % IEC 61400-1 NTM.
%                                      % Class B, Iref = 0.14.

sigV = 1.;   % Unit variance.
TI = sigV/Vmag;

%Lu    = (180/length); % Based on surface roughness length of
%                      % about 0.01 m.  0.0001: 235, 0.001: 200,
%                      % 0.01: 175, 0.1: 150.  Intended to
%                      % represent a turbine deep in a plant.

Lu = 2000;

nnf = 2^9;  % Number of analysis frequencies, 1/4 the total frequencies.
df = 0.001;

%-------------------------------------------------------------------

Pin = assemblePin (s);
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);
[q0,P,Ts0B,TB0g] =                                                     ...
      undeformedPosition (Pin,yaw,s.nacelle.delta,0,s.driveshaft.phi, ...
                          betas,0,idofs,idofm,inods,inodm);
Ndj = size(q0,1);

%-------------------------------------------------------------------
% Specify the wind shear as a function of height.
% Elevation, vector pairs.
Nws = 50;
z1 = P(idofs(2)+3) - 10;
z2 = P(idofs(4)+3) + 120;
zref = P(idofs(4)+3) - z1;
h0 = 0.01;
dz = (z2 - z1 - h0)/(Nws - 1);
zs = h0 + dz*[0:Nws-1].';
WS = (Vinf.').*log(zs/h0)/log(zref/h0);
%-------------------------------------------------------------------

ppx = pchip (zs+z1,WS(:,1));
ppy = pchip (zs+z1,WS(:,2));
ppz = pchip (zs+z1,WS(:,3));

Nf = 4*nnf;
dt = 1/(4*nnf*df);

[b1,b2,b3] = MBCindices_Ndj (Ndj,idofs);
[TpsiB_Ndj,TBpsi_Ndj] = MBC (Ndj,b1,b2,b3,psi0);
q = TpsiB_Ndj*ypsi(1:Ndj);
dqdt = TpsiB_Ndj*ypsi(Ndj+[1:Ndj]);

QSflag = 2;
Sij = buildSij (QSflag,MBCflag,s,a,nnf,df,Vinf,TI,Lu,q,dqdt,P);

Npsi = 2^11;
psis = (2*pi/Npsi)*[0:Npsi-1].';

[Vavg,vg] = periodicWind (s,a,q,P,psis,ppx,ppy,ppz);

% vg is reported as a function of the azimuth angle psi.  I need
% this in terms of time, psi/W.
Omega = ypsi(Ndj+idofs(4)+6);

ts = psis/Omega;
dt = ts(2) - ts(1);
P1 = Omega/(2*pi);    % 1P = df = Omega/(2*pi) = 1/(Npsi*dt).

% Also the MBC transform.
Nel = size(vg,2)/3;

if (MBCflag == 1)
   TBpsi = zeros (Npsi,9);
   third = 1/3;
   cp  = cos (psis);
   sp  = sin (psis);
   cp2 = cos (psis + 2*pi/3);
   sp2 = sin (psis + 2*pi/3);
   cp3 = cos (psis + 4*pi/3);
   sp3 = sin (psis + 4*pi/3);
   TBpsi(:,1) = third;
   TBpsi(:,2) = 2*third*cp;
   TBpsi(:,3) = 2*third*sp;
   TBpsi(:,4) = third;
   TBpsi(:,5) = 2*third*cp2;
   TBpsi(:,6) = 2*third*sp2;
   TBpsi(:,7) = third;
   TBpsi(:,8) = 2*third*cp3;
   TBpsi(:,9) = 2*third*sp3;
   ibl = [[1:Nel];Nel+[1:Nel];2*Nel+[1:Nel]];
   for iel = 1:Nel  % Note, really iel = icol, indexing the component.
      ir3 = 3*(iel-1);
      dat = vg(:,ibl(:,iel)) + Vavg(ibl(:,iel)).';  % Need to include mean.
      vv = zeros(Npsi,3);
      vv(:,1) = TBpsi(:,1).*dat(:,1) + TBpsi(:,4).*dat(:,2) ...
              + TBpsi(:,7).*dat(:,3);
      vv(:,2) = TBpsi(:,2).*dat(:,1) + TBpsi(:,5).*dat(:,2) ...
              + TBpsi(:,8).*dat(:,3);
      vv(:,3) = TBpsi(:,3).*dat(:,1) + TBpsi(:,6).*dat(:,2) ...
              + TBpsi(:,9).*dat(:,3);
      Vavg(ibl(:,iel)) = mean(vv).';
      vg(:,ibl(1,iel)) = vv(:,1) - Vavg(ibl(1,iel));
      vg(:,ibl(2,iel)) = vv(:,2) - Vavg(ibl(2,iel));
      vg(:,ibl(3,iel)) = vv(:,3) - Vavg(ibl(3,iel));
   end
end

vgf = fft(vg)/Npsi;
Nv = size(vg,2);

Sper = zeros (Npsi,Nv^2);
for iv = 1:Nv

   icol = Nv*(iv-1);
   Sper(:,icol+[1:Nv]) = vgf.*conj(vgf(:,iv))/df;

end

% Store the spectra according to the actual frequencies, not multiples
% of the rotor frequency.
ipf = round(([0:Npsi/2-1].')*P1/df);
jpf = ipf<(Nf/2-1);
ipfr = ipf(jpf) + 1;

Sijp = sparse(Nf,Nv^2);
Sijp(ipfr,:) = Sper(jpf,:);

temp = vgf;
vgf = sparse(Nf,Nv);
vgf(ipfr,:) = temp(jpf,:);

jpsi = [-Npsi/2:-1].';
ipf = round(jpsi*P1/df);
jpf = ipf>(-Nf/2);
ipfr = ipf(jpf);
Sijp(Nf+ipfr+1,:) = Sper(Npsi+jpsi(jpf)+1,:);
vgf(Nf+ipfr+1,:) = temp(Npsi+jpsi(jpf)+1,:);

cvg = vgf;

save ('-float-binary','Sij.bin','Sij');
save ('-float-binary','Sijp.bin','Sijp');
save ('-float-binary','cvg.bin','cvg');
save ('-float-binary','Vavg.bin','Vavg');


