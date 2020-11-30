
clear;

pkg load statistics;

nm = 'DTU10MW';
eval(['[s,a] = STASTurbine_'  nm ' ();']);
eval(["load 'cpct_" nm ".txt';"]);
eval(["cpct = cpct_" nm ";"]);

[length,time,mass,current,voltage,              ...
 velocity,force,power,stress,density,viscosity, ...
 stiffness,damping,resistance,inductance,       ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
Ndj = Ndof + 6;

outnm = '_ConstP';

% Load subsystem matrices.
inpnm = '_P060_V100';
eval(["load 'AT"   inpnm ".bin';"]);   % Turbines, discrete-time
eval(["load 'BT"   inpnm ".bin';"]);
eval(["load 'Cq"   inpnm ".bin';"]);
eval(["load 'Dq"   inpnm ".bin';"]);
eval(["load 'AG"   inpnm ".bin';"]);   % Observers, discrete-time
eval(["load 'BG"   inpnm ".bin';"]);
eval(["load 'GG"   inpnm ".bin';"]);
eval(["load 'KI"   inpnm ".bin';"]);
eval(["load 'CG"   inpnm ".bin';"]);
eval(["load 'DG"   inpnm ".bin';"]);
eval(["load 'xpsi" inpnm ".txt';"]);
eval(["load 'ypsi" inpnm ".txt';"]);
eval(["xpsi = xpsi" inpnm ";"]);
eval(["ypsi = ypsi" inpnm ";"]);

inpnm = '_TCRWP_P060';
eval(["load 'Bgq"  inpnm ".bin';"]);   % Quasi-steady grid.
eval(["load 'gret" inpnm ".bin';"]);
eval(["load 'xgf" inpnm ".bin';"]);

xgf0 = xgf;  % Reserve xgf for the perturbations.

Vref = 10                              / velocity;
Dia  = 2*89.15                         / length;
we   = 50*2*pi                         * time;    % [rad/s], grid frequency
vpcc = [400;0].'                       / voltage;

% Wind cascade.
dtvC = 1                               / time;
NvC  = 2^8;  % That'll be a bit over 4 minutes.

% Fatigue.
%dtF = 10                               / time;

% Cluster probabilities.
gam  = 0.01;                  % Binary white noise amplitude.
sig0 = 2.0;                   % Gaussian stdev.
bet0 = 0.2;                   % Background probability.
dmu  = 0.1;
mbnd = [-0.5*Vref dmu 0.5*Vref];

% Plant power and other references.
P0 = 6.e6                              / power;
Ppcc = 32*P0;
R = 0.5*Dia;
W0 = ypsi(Ndj+idofs(4)+6);
pitch0 = -ypsi(idofs(6)+4);
[cp,ct,dcp,dct] = cpvwb (cpct,R,Vref,W0,pitch0);
Area = pi*(R^2);
dens = 1.225                           / density;
FT0 = ct*0.5*dens*Area*(Vref^2);
uc = [Ppcc, Ppcc, Vref, P0, FT0, W0, 0, 31, 31].';
xc0 = [0;0;zeros(3*9,1)];
xc0(3+[0:3:3*(9-1)]) = uc;



global betaglobal;
betaglobal = pitch0;



[cp,ct,dcp,dct] = cpvwb (cpct,R,Vref,W0,0);
Pa0 = cp*0.5*dens*Area*(Vref^3);

a1table = interp1 ([0, 0.25,  0.50,  0.75,   1], ...
                   [1, 0.975, 0.95,  0.925,  0.90],'pchip','pp');
a2table = interp1 ([0, 0.25,  0.50,  0.75,   1], ...
                   [1, 0.75,  0.50,  0.25,   0.00],'pchip','pp');
atables.a1 = a1table;
atables.a2 = a2table;

p = [1.0;
     0.4                                * time; ...
     6.0                                / velocity; ...
     0.03*(2*pi)                        * time; ...
    -1.e6                               / power; ...
     1.e6                               / power; ...
     0.10; ...
     Area; ...
     dens; ...
     0.08*(2*pi)                        * time; ...
     0.24*(2*pi)                        * time; ...
     0.04; ...
     0.40];

% Simulation.
dt = 0.05                               / time;
Nt = 72000;
ts = dt*[0:Nt-1].';

dtsave = 0.25                           / time;
Ntsave = Nt/round(dtsave/dt);

NtvC = Nt/round(dtvC/dt);
Nf = 2*NvC;
df = 1/(Nf*dtvC);
freqs = df*[[0:Nf/2-1] [-Nf/2:-1]].';
tvCs = dtvC*[0:NtvC-1].';

%NtF = Nt/round(dtF/dt);
%tFs = dtF*[0:NtF-1].';

% Load turbine positions.
load 'TurbinePosition.dat';  % x/D, y/D.
npos = TurbinePosition*Dia;
Nturb = size(npos,1);

% Load wind files.  W: Vmag, thV for each turbine.
% Truncate to 3 hours, eliminating the earliest timesteps.  Valentin
% said there may be a "burn-in" period at the beginning.
dtW = 5                                / time;
NtW = 1801;
tWs = dtW*[0:NtW-1];
load 'WindTS.bin';
NW = size(W,1);
W = W(NW-NtW+[1:NtW],:)                / velocity;

%load 'Sij_V10_I183_Lu180_W911.bin';
%load 'Sijp_V10_I183_Lu180_W911.bin';
%load 'dSdV_V10_I183_Lu180_W911.bin';  % For special method of scaling with W.
%TInom = 0.183;
%Wnom = 0.91146;  % Rotor speed used when generating this particular turbulence file.

Nxt = size(AT,1);
Nyt = size(Cq,1);
Nut = size(BT,2);
Nxo = size(AG,1);
Nyo = size(CG,1);
Nuo = size(BG,2);
Nso = size(GG,2);
Nxc = 29;
Nyc = 1;
Nuc = 9;
Nxi = Nso;

iusens = [4, 2, 3];        % Pe, beta, yaw.
iysens = [1, 5, 6, 7, 8];  % Omega, vnacx,y, V, thV.

% The full grid states.  Bgq is based on the reduced "gret" DOFs.
% Inputs to Bgq: we(1), vPCCd,q(2), itd,q(2*Ngnod)
Ngcab = 32;
Ngnod = 33;  % Collector grid and substation.
Nxgf = 2*Ngcab + 4*Ngnod + 14;  % Cable currents, bus voltages, bus shunt currents,
                                % and the export system states.

% Determine clusters.
[dcl,icl] = cluster (npos.',9);

% Initialize the cascade with the wind speed at t = 0.
vcas = W(1,1:2:2*Nturb-1).*ones(NvC,Nturb) - Vref;
[AvC,BvC] = Vcascade (NvC);

% Initialize the cluster probabilities.
sig = sig0*ones(Nturb,1);
bet = bet0*ones(Nturb,1);
mus = [mbnd(1):mbnd(2):mbnd(3)].';
Nmu = size(mus,1);
phi = (1/Nmu)*ones(Nmu,Nturb);  % Uniform prior.
Ncells = Nmu;
ics = [1:Ncells].';
mc = xFromCell (ics,mbnd);
mc1 = mc + gam*dt;
[ics,wcs] = nearestCells (mc1,mbnd);
Lam = assignPhi (ics,0.5*wcs,Ncells);
mc1 = mc - gam*dt;
[ics,wcs] = nearestCells (mc1,mbnd);
L2 = assignPhi (ics,0.5*wcs,Ncells);
Lam = Lam + L2;  % Transition probability matrix.

% Step in time; use the smallest timestep, that of the observer
% and turbine models.
del = sqrt(eps);

xpc = zeros(Ntsave,2*Nturb);
yts = zeros(Ntsave,11*Nturb);
yos = zeros(Ntsave,11*Nturb);
Phs = zeros(Ntsave,Nturb);
ucs = zeros(Ntsave,Nuc*Nturb);
Vs  = zeros(Ntsave,Nturb);

xt = zeros(Nxt*Nturb,1);
xo = zeros(Nxo*Nturb,1);
xc = zeros(Nxc*Nturb,1);
for iturb = 1:Nturb
   xc(Nxc*(iturb-1)+[1:Nxc]) = xc0;
end
xi = zeros(Nxi*Nturb,1);

yt = zeros(Nyt*Nturb,1);
yo = zeros(Nyo*Nturb,1);
Ph = P0*ones(Nturb,1);

uc = zeros(Nuc*Nturb,1);

Dam = zeros(Nturb,1);
a1sum = zeros(Nturb,1);
a2sum = zeros(Nturb,1);

rand('state',1234);
Dam = rand(Nturb,1);



return



%Dam(1) = 1;
ind = [1:Nturb].';
alf1 = zeros(Nturb,1);
alf2 = zeros(Nturb,1);
for iturb = 1:Nturb
   [alf1(iturb),da1dD] = gains1 (Dam(iturb),a1table);
   [alf2(iturb),da2dD] = gains1 (Dam(iturb),a2table);
end
for iturb = 1:Nturb
   a1sum(iturb) = sum(alf1(ind ~= iturb));
   a2sum(iturb) = sum(alf2(ind ~= iturb));
end

iPhPCC = 1 + [0:9:9*(Nturb-1)].';
iPmPCC = 2 + [0:9:9*(Nturb-1)].';
imuV   = 3 + [0:9:9*(Nturb-1)].';
iPa    = 4 + [0:9:9*(Nturb-1)].';
iFT    = 5 + [0:9:9*(Nturb-1)].';
iW     = 6 + [0:9:9*(Nturb-1)].';
iDam   = 7 + [0:9:9*(Nturb-1)].';
ia1sum = 8 + [0:9:9*(Nturb-1)].';
ia2sum = 9 + [0:9:9*(Nturb-1)].';

iyW   = 1  + [0:Nyt:Nyt*(Nturb-1)].';
iybet = 2  + [0:Nyt:Nyt*(Nturb-1)].';
iyyaw = 3  + [0:Nyt:Nyt*(Nturb-1)].';
iyPe  = 4  + [0:Nyt:Nyt*(Nturb-1)].';
iyvnx = 5  + [0:Nyt:Nyt*(Nturb-1)].';
iyvny = 6  + [0:Nyt:Nyt*(Nturb-1)].';
iyV   = 7  + [0:Nyt:Nyt*(Nturb-1)].';
iythV = 8  + [0:Nyt:Nyt*(Nturb-1)].';
iyFT  = 10 + [0:Nyt:Nyt*(Nturb-1)].';
iyPa  = 11 + [0:Nyt:Nyt*(Nturb-1)].';

iistd = 60 + [0:Nyt:Nyt*(Nturb-1)].';
iistq = 61 + [0:Nyt:Nyt*(Nturb-1)].';
iist = sort([iistd;iistq],'ascend');

ivstd = 2*Ngcab + [1:2:2*(Ngnod-1)].';
ivstq = 2*Ngcab + [2:2:2*Ngnod].';
ivst = sort([ivstd;ivstq],'ascend');

eval(["fidp = fopen ('FTP" outnm ".txt','w');"]);
eval(["fidv = fopen ('VV" outnm ".txt','w');"]);

itvC = 1;
%itF = 1;
itsave = 0;
Cflg = false;
%Fflg = false;
Sflg = false;
for it = 1:Nt

   t = ts(it);
   if (abs(mod(t,dtvC)) < del) || (abs(mod(t,dtvC)-1) < del)
      Cflg = true;
      itvC = itvC + 1;
      ict = Nturb*(itvC-1);
printf('%10.3f\n',t);
fflush(stdout);
   else
      Cflg = false;
   end

%   if (abs(mod(t,dtF)) < del) || (abs(mod(t,dtF)-1) < del)
%      Fflg = true;
%      itF = itF + 1;
%   else
%      Fflg = false;
%   end

   if (abs(mod(t,dtsave)) < del) || (abs(mod(t,dtsave)-1) < del)
      Sflg = true;
      itsave = itsave + 1;
      fprintf(fidp,'%+5.6e',itsave*dtsave);
      fprintf(fidv,'%+5.6e',itsave*dtsave);
   else
      Sflg = false;
   end

   PhPCC = Ppcc;

   % Update the wind inputs.
   V = interp1 (tWs,W,t);
   dV = zeros(2*Nturb,1);
   dV(1:2:2*Nturb-1) = V(1:2:2*Nturb-1) - Vref;
   dV(2:2:2*Nturb)   = V(2:2:2*Nturb);

   if (Sflg)
      Vs(itsave,:) = V(1:2:2*Nturb-1);
   end

   % Update the grid model.  we, vpccd,q, turbine currents (d,q).
   ug = [0;0;0;0.5*yt(iist);0;0];  % 0.5: 33 kV turbine model, 66 kV grid.
   xg = Bgq*ug;
   xgf = sparse(Nxgf,1);
   xgf(gret) = xg;
   iref = 2*Ngcab + 4*Ngnod;
   PmPCC = Ppcc                                                    ...
         + xgf0(iref+11)*xgf(iref+13) + xgf0(iref+12)*xgf(iref+14) ...
         + xgf(iref+11)*xgf0(iref+13) + xgf(iref+12)*xgf0(iref+14);
   vst = xgf(ivst);

   for iturb = 1:Nturb

      icxt = Nxt*(iturb-1);
      icxo = Nxo*(iturb-1);
      icxc = Nxc*(iturb-1);
      icxi = Nxi*(iturb-1);
      icuc = Nuc*(iturb-1);
      icyt = Nyt*(iturb-1);
      icyo = Nyo*(iturb-1);
      ic11 = 11*(iturb-1);
      ic2 = 2*(iturb-1);

      % Update the dynamic turbine models.
%      ut = [Ph(iturb)-P0;dV(ic2+[1:2]);0;vst(ic2+[1:2])];
ut = [Ph(iturb)-P0;dV(ic2+[1:2]);0;0;0];  % Ignore terminal voltage, numerically unstable.
      xt(icxt+[1:Nxt]) = AT*xt(icxt+[1:Nxt]) + BT*ut;
      yt(icyt+[1:Nyt]) = Cq*xt(icxt+[1:Nxt]) + Dq*ut;

      if (Sflg)
         yts(itsave,ic11+[1:11]) = yt(icyt+[1:11]);
         fprintf(fidp,' %+5.6e %+5.6e',FT0+yt(icyt+10),P0+yt(icyt+4));
      end

      % Update the observer models.
      uo = yt(icyt+iusens);
      xo(icxo+[1:Nxo]) = AG*xo(icxo+[1:Nxo]) + BG*uo ...
                       + GG*(yt(icyt+iysens) + xi(icxi+[1:Nxi]));
      yo(icyo+[1:Nyo]) = CG*xo(icxo+[1:Nxo]) + DG*uo;
      xi(icxi+[1:Nxi]) = xi(icxi+[1:Nxi]) ...
                       + KI*(yt(icyt+iysens) - yo(icyo+iysens));

      if (Sflg)
         yos(itsave,ic11+[1:11]) = yo(icyo+[1:11]);
      end

      % Update the cluster wind speed estimates and cascade.
      if (Cflg)

         vobs = yo(iyV);
         phi0 = phi(:,iturb);
         Vdat = vobs(icl(:,iturb));  % The cluster.
         [dmuV,phi(:,iturb)] = ...
            clusdt (phi0,Vdat,Lam,mus,sig(iturb),bet(iturb));
         vcas(:,iturb) = AvC*vcas(:,iturb) ...
                       + BvC*(vobs(iturb) - dmuV);
         uc(imuV(iturb)) = dmuV + Vref;

         % Scale sig to a reasonable value with respect to the present
         % wind speed.
         sig(iturb) = 0.2*uc(imuV(iturb));

      end

      if (Sflg)
         fprintf(fidv,' %+5.6e %+5.6e',V(2*(iturb-1)+1),uc(imuV(iturb)));
      end

      % Update the estimate of fatigue.  [Leave this for later.]
%      if (Fflg)

         % Get the first two spectral moments of the turbulence.

         % Estimate VK spectra parameters based on the spectral moments.

         % Scale the reference spectrum to the current W, I, Lu.

%      end


      % Update the plant controller.
      uc(iPhPCC(iturb)) = PhPCC;
      uc(iPmPCC(iturb)) = PmPCC;
      % (muV is input above.)
      uc(iPa(iturb))    = Pa0 + yt(iyPa(iturb));
      uc(iFT(iturb))    = FT0 + yt(iyFT(iturb));
      uc(iW(iturb))     = W0  + yt(iyW(iturb));
      uc(iDam(iturb))   = Dam(iturb);
      uc(ia1sum(iturb)) = a1sum(iturb);
      uc(ia2sum(iturb)) = a2sum(iturb);
%{
      [dxdt,Ph(iturb),~,~,~,~] = plantControl (xc(icxc+[1:Nxc]),uc(icuc+[1:Nuc]), ...
                                               p,atables,cpct,0);
      x1 = xc(icxc+[1:Nxc]) + 0.5*dt*dxdt;
      [dxdt,Ph(iturb),~,~,~,~] = plantControl (x1,uc(icuc+[1:Nuc]), ...
                                               p,atables,cpct,0);
      xc(icxc+[1:Nxc]) = xc(icxc+[1:Nxc]) + dt*dxdt;
%}
Ph(iturb) = P0;

      if (Sflg)
         xpc(itsave,ic2+[1:2]) = xc(icxc+[1:2]);
         Phs(itsave,iturb) = Ph(iturb);
         ucs(itsave,icuc+[1:Nuc]) = uc(icuc+[1:Nuc]);
      end

   end  % Turbine

   if (Sflg)
      fprintf(fidp,'\n');
      fprintf(fidv,'\n');
   end

end  % Time

clrz = [0.0, 0.0, 0.0; ...
        0.6, 0.6, 0.6; ...
        0.8, 0.0, 0.0; ...
        0.4, 0.0, 0.0; ...
        0.0, 0.8, 0.0; ...
        0.0, 0.4, 0.0; ...
        0.0, 0.0, 0.8; ...
        0.0, 0.0, 0.4];

iturb = 1;

%{
% Control block interface variables.
figure(1);
clf;
hold on;
plot(ts,Phs(:,iturb),'color',clrz(6,:),'linewidth',2);          % Phat  (dark green)
plot(ts,ucs(:,imuV(iturb)),'color',clrz(1,:),'linewidth',2);    % muV   (black)
plot(ts,ucs(:,iPmPCC(iturb))/100,'color',clrz(5,:),'linewidth',2);  % PmPCC (light green)
plot(ts,ucs(:,iFT(iturb)),'color',clrz(2,:),'linewidth',2);     % FT    (gray)
plot(ts,ucs(:,iPa(iturb)),'color',clrz(3,:),'linewidth',2);     % Pa    (light red)
plot(ts,ucs(:,iW(iturb)),'color',clrz(7,:),'linewidth',2);      % W     (blue)
hold off;

i2a = [1:2:2*Nturb-1].';
i2b = [2:2:2*Nturb].';
figure(2);
clf;
hold on;
plot(ts,Phs(:,iturb),'color',clrz(6,:),'linewidth',2);          % Phat  (dark green)
plot(ts,xpc(:,i2a(iturb)),'color',clrz(3,:),'linewidth',2);     % intP  (light red)
plot(ts,xpc(:,i2b(iturb)),'color',clrz(4,:),'linewidth',2);     % intT  (dark red)
hold off;

% Turbine.
figure(3);
clf;
hold on;
plot(ts,Phs(:,iturb),'color',clrz(6,:),'linewidth',2);          % Phat  (dark green)
plot(ts,Vs(:,iturb),'color',clrz(1,:),'linewidth',2);           % V     (black)
plot(ts,P0+yts(:,11*(iturb-1)+iyPe(1)),'color',clrz(5,:),'linewidth',2); % Pe (lgt grn)
plot(ts,ucs(:,iW(iturb)),'color',clrz(7,:),'linewidth',2);      % W     (blue)
hold off;

% Observer wind output versus actual.
figure(4);
clf;
hold on;
plot(ts,Vref+yos(:,11*(iturb-1)+iyV(1)),'color',clrz(3,:),'linewidth',2);  % Vobs (light red)
plot(ts,Vref+yts(:,11*(iturb-1)+iyV(1)),'color',clrz(4,:),'linewidth',2);  % V    (dark red)
hold off;
%}

fclose('all');