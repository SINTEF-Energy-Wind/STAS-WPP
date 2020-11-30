% y vector:        u vector:
% (controls)       muV      1
% Phat     1       Dam      2
% Pe       2       asum     3
% beta     3       PhPCC    4
% yaw      4       PmPCC    5
% (sensors)        V        6  % Rotor-avg V at turbine.
% W        5
% vnacx,y 6-7
% V,thV   8-9
% (estimates)
% Vst     10
% FTst    11
% Past    12
% (grid)
% ist   13,14

clear;

pkg load statistics;

nm = 'DTU10MW';
eval(['[s,a] = STASTurbine_'  nm ' ();']);

[length,time,mass,current,voltage,        ...
 velocity,force,power,stress,ndens,nvisc, ...
 stiffness,damping,resistance,inductance, ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
Ndj = Ndof + 6;

Vref = 10                              / velocity;
Dia  = 2*89.15                         / length;
we   = 50*2*pi                         * time;    % [rad/s], grid frequency
vpcc = [400;0].'                       / voltage;

% Wind cascade.
dtvC = 1                               / time;
NvC  = 2^8;  % That'll be about 4 minutes.

% Cluster probabilities.
gam  = 0.01;                  % Binary white noise amplitude.
sig0 = 2.0;                   % Gaussian stdev.
bet0 = 0.2;                   % Background probability.
dmu  = 0.1;
%mbnd = [-0.5*Vref dmu 0.5*Vref];
mbnd = [-Vref dmu Vref];

% Load turbine positions.
load 'TurbinePosition.dat';  % x/D, y/D.
npos = TurbinePosition*Dia;
Nturb = size(npos,1);

a1table = interp1 ([0, 0.25,  0.50,  0.75,   1], ...
                   [1, 0.975, 0.95,  0.925,  0.90],'pchip','pp');
a2table = interp1 ([0, 0.25,  0.50,  0.75,   1], ...
                   [1, 0.75,  0.50,  0.25,   0.00],'pchip','pp');
atables.a1 = a1table;
atables.a2 = a2table;

% Damage and production factors, not yet updated dynamically.
Dam = zeros(Nturb,1);
a1sum = zeros(Nturb,1);
a2sum = zeros(Nturb,1);

%iturb = 24;
%Dam(iturb) = 1;

rand('state',1234);
%Dam = rand(Nturb,1);
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

% Simulation.
dt = 0.05                              / time;
Nt = 72000;
ts = dt*[0:Nt-1].';

NtvC = Nt/round(dtvC/dt);
Nf = 2*NvC;
df = 1/(Nf*dtvC);
freqs = df*[[0:Nf/2-1] [-Nf/2:-1]].';
tvCs = dtvC*[0:NtvC-1].';

% Load wind files.  W: Vmag, thV for each turbine.
% Truncate to 3 hours, eliminating the earliest timesteps.  Valentin
% said there may be a "burn-in" period at the beginning.
dtW = 5                                / time;
NtW = 1801;
tWs = dtW*[0:NtW-1];
load 'WindTS.bin';
NW = size(W,1);
W = W(NW-NtW+[1:NtW],:)                / velocity;

% Load subsystem matrices.
inpnm = '_P060_V100';
eval(["load 'ACL"   inpnm ".bin';"]);
eval(["load 'BCL"   inpnm ".bin';"]);
eval(["load 'CCL"   inpnm ".bin';"]);
eval(["load 'DCL"   inpnm ".bin';"]);
eval(["load 'xpsi" inpnm ".txt';"]);
eval(["load 'ypsi" inpnm ".txt';"]);
eval(["xpsi = xpsi" inpnm ";"]);
eval(["ypsi = ypsi" inpnm ";"]);


inpnm = '_TCRWP_P060';
eval(["load 'Bgq"  inpnm ".bin';"]);   % Quasi-steady grid.
eval(["load 'gret" inpnm ".bin';"]);
eval(["load 'xgf" inpnm ".bin';"]);

xgf0 = xgf;  % Reserve xgf for the perturbations.

%load 'Sij_V10_I183_Lu180_W911.bin';
%load 'Sijp_V10_I183_Lu180_W911.bin';
%load 'dSdV_V10_I183_Lu180_W911.bin';  % For special method of scaling with W.
%TInom = 0.183;
%Wnom = 0.91146;  % Rotor speed used when generating this particular turbulence file.

Nx = size(ACL,1);
Nu = size(BCL,2);
Ny = size(CCL,1);

Nxt = 44;
Nxo = 25;
Nxc = 29;

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
%xpc = zeros(Nt,2*Nturb);
%ys = zeros(Nt,Ny*Nturb);
%us = zeros(Nt,Nu*Nturb);

x = zeros(Nx*Nturb,1);
y = zeros(Ny*Nturb,1);
u = zeros(Nu*Nturb,1);

i2a = 1 + [0:2:2*(Nturb-1)].';
i2b = 2 + [0:2:2*(Nturb-1)].';

imuV = 1 + [0:7:7*(Nturb-1)].';
iDam = 2 + [0:7:7*(Nturb-1)].';
ia1s = 3 + [0:7:7*(Nturb-1)].';
ia2s = 4 + [0:7:7*(Nturb-1)].';
iPh  = 5 + [0:7:7*(Nturb-1)].';
iPm  = 6 + [0:7:7*(Nturb-1)].';
iV   = 7 + [0:7:7*(Nturb-1)].';

iPhat = 1  + [0:14:14*(Nturb-1)].';
iPe   = 2  + [0:14:14*(Nturb-1)].';
ibet  = 3  + [0:14:14*(Nturb-1)].';
iyaw  = 4  + [0:14:14*(Nturb-1)].';
iW    = 5  + [0:14:14*(Nturb-1)].';
ivnx  = 6  + [0:14:14*(Nturb-1)].';
ivny  = 7  + [0:14:14*(Nturb-1)].';
iVV   = 8  + [0:14:14*(Nturb-1)].';
ithV  = 9  + [0:14:14*(Nturb-1)].';
iVst  = 10 + [0:14:14*(Nturb-1)].';
iFTst = 11 + [0:14:14*(Nturb-1)].';
iPast = 12 + [0:14:14*(Nturb-1)].';
iistd = 13 + [0:14:14*(Nturb-1)].';
iistq = 14 + [0:14:14*(Nturb-1)].';
iist = sort([iistd;iistq],'ascend');

fidv = fopen ('tcas.txt','w');

itvC = 1;
Cflg = false;
for it = 1:Nt

   t = ts(it);
   if (abs(mod(t,dtvC)) < del)
      Cflg = true;
      itvC = itvC + 1;
      ict = Nturb*(itvC-1);
printf('%12.3f\n',t);
fflush(stdout);
   else
      Cflg = false;
   end

   PhPCC = 0;

   % Update the wind inputs.
   V = interp1 (tWs,W,t);



%{
if (it < 12000)
V = interp1 (tWs,W,t);
else
V = interp1 (tWs,W,t);
V = V - 5;
end
V(2*11+1) = 0;
V(2*18+1) = 0;
%}



   dV = zeros(2*Nturb,1);
   dV(1:2:2*Nturb-1) = V(1:2:2*Nturb-1) - Vref;
   dV(2:2:2*Nturb)   = V(2:2:2*Nturb);
   u(iV)   = dV(1:2:2*Nturb-1);

   % Update the grid model.  we, vpccd,q, turbine currents (d,q).
   ug = [0;0;0;0.5*y(iist);0;0];  % 0.5: 33 kV turbine model, 66 kV grid.
   xg = Bgq*ug;
   xgf = sparse(Nxgf,1);
   xgf(gret) = xg;
   iref = 2*Ngcab + 4*Ngnod;
   u(iPm) = xgf0(iref+11)*xgf(iref+13) + xgf0(iref+12)*xgf(iref+14) ...
          + xgf(iref+11)*xgf0(iref+13) + xgf(iref+12)*xgf0(iref+14);

   u(iDam) = Dam;  % Constant for now, until the spectral part is working.
   u(ia1s) = a1sum;
   u(ia2s) = a2sum;
   u(iPh)  = PhPCC;

   fprintf (fidv,'%+5.6e',t);

   for iturb = 1:Nturb

      ic2 = 2*(iturb-1);
      ic6 = 6*(iturb-1);
      icx = Nx*(iturb-1);
      icu = Nu*(iturb-1);
      icy = Ny*(iturb-1);

      % Update the dynamic turbine models.
      % Inputs: muV, Dam, asum, PhPCC, PmPCC, V
      x(icx+[1:Nx]) = ACL*x(icx+[1:Nx]) + BCL*u(icu+[1:Nu]);
      y(icy+[1:Ny]) = CCL*x(icx+[1:Nx]) + DCL*u(icu+[1:Nu]);

      % Update the cluster wind speed estimates and cascade.
      if (Cflg)

         vobs = y(iVst);
         phi0 = phi(:,iturb);
         Vdat = vobs(icl(:,iturb));  % The cluster.
         [u(imuV(iturb)),phi(:,iturb)] = ...
            clusdt (phi0,Vdat,Lam,mus,sig(iturb),bet(iturb));
         vcas(:,iturb) = AvC*vcas(:,iturb) ...
                       + BvC*(vobs(iturb) - u(imuV(iturb)));

         % Scale sig to a reasonable value with respect to the present
         % wind speed.
         sig(iturb) = 0.2*(u(imuV(iturb)) + Vref);

      end

%      xpc(it,ic2+[1:2]) = x(icx+Nxt+Nxo+[1:2]);
%      ys(it,icy+[1:Ny]) = y(icy+[1:Ny]);
%      us(it,icu+[1:Nu]) = u(icu+[1:Nu]);

      fprintf (fidv,' %+5.6e %+5.6e',V(ic2+1),u(imuV(iturb))+Vref);

   end % iturb

   fprintf (fidv,'\n');

end % it

clrz = [0.0, 0.0, 0.0; ...
        0.6, 0.6, 0.6; ...
        0.8, 0.0, 0.0; ...
        0.4, 0.0, 0.0; ...
        0.0, 0.8, 0.0; ...
        0.0, 0.4, 0.0; ...
        0.0, 0.0, 0.8; ...
        0.0, 0.0, 0.4];

%{

iturb = 2;

% Control block interface variables.
figure(1);
clf;
hold on;
plot(ts,ys(:,iPhat(iturb)),'color',clrz(6,:),'linewidth',2);
plot(ts,us(:,imuV(iturb)),'color',clrz(1,:),'linewidth',2);
plot(ts,us(:,iPm(iturb))/10,'color',clrz(5,:),'linewidth',2);
plot(ts,ys(:,iFTst(iturb)),'color',clrz(2,:),'linewidth',2);
plot(ts,ys(:,iPast(iturb)),'color',clrz(3,:),'linewidth',2);
plot(ts,ys(:,iW(iturb)),'color',clrz(7,:),'linewidth',2);
hold off;

figure(2);
clf;
hold on;
plot(ts,ys(:,iPhat(iturb)),'color',clrz(6,:),'linewidth',2);
plot(ts,xpc(:,i2a(iturb)),'color',clrz(3,:),'linewidth',2);
plot(ts,xpc(:,i2b(iturb)),'color',clrz(4,:),'linewidth',2);
hold off;

% Turbine.
figure(3);
clf;
hold on;
plot(ts,ys(:,iPhat(iturb)),'color',clrz(6,:),'linewidth',2);
plot(ts,us(:,iV(iturb)),'color',clrz(1,:),'linewidth',2);
plot(ts,ys(:,iPe(iturb)),'color',clrz(5,:),'linewidth',2);
plot(ts,ys(:,iW(iturb)),'color',clrz(7,:),'linewidth',2);
hold off;

% Observer wind output versus actual.
figure(4);
clf;
hold on;
plot(ts,ys(:,iVst(iturb)),'color',clrz(3,:),'linewidth',2);
plot(ts,us(:,iV(iturb)),'color',clrz(4,:),'linewidth',2);
hold off;

%}

fclose('all');