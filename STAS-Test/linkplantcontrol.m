
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

% Load subsystem matrices.
inpnm = '_P060_V100';
eval(["load 'AT"   inpnm ".bin';"]);
eval(["load 'BT"   inpnm ".bin';"]);
eval(["load 'Cq"   inpnm ".bin';"]);
eval(["load 'Dq"   inpnm ".bin';"]);
eval(["load 'AG"   inpnm ".bin';"]);
eval(["load 'BG"   inpnm ".bin';"]);
eval(["load 'CG"   inpnm ".bin';"]);
eval(["load 'DG"   inpnm ".bin';"]);
eval(["load 'GG"   inpnm ".bin';"]);
eval(["load 'xt"   inpnm ".bin';"]);
eval(["load 'ut"   inpnm ".bin';"]);
eval(["load 'yt"   inpnm ".bin';"]);
eval(["load 'xo"   inpnm ".bin';"]);
eval(["load 'uo"   inpnm ".bin';"]);
eval(["load 'yo"   inpnm ".bin';"]);
eval(["load 'xc"   inpnm ".bin';"]);
eval(["load 'uc"   inpnm ".bin';"]);
eval(["load 'Phat" inpnm ".bin';"]);
eval(["load 'xpsi" inpnm ".txt';"]);
eval(["load 'ypsi" inpnm ".txt';"]);
eval(["xpsi = xpsi" inpnm ";"]);
eval(["ypsi = ypsi" inpnm ";"]);


Nxt = size(AT,1);
Nut = size(BT,2);
Nyt = size(Cq,1);

Nxo = size(AG,1);
Nuo = size(BG,2);
Nyo = size(CG,1);

% Mean values.
Vinf = 10                               / velocity;
P0 = Phat;
bet0 = -ypsi(idofs(6)+4) + yt(2); 

R = 89.15                               / length;
W0 = uc(6);
[cp,ct,dcp,dct] = cpvwb (cpct,R,Vinf,W0,bet0);
Area = pi*(R^2);
dens = 1.225                            / density;
FT0 = uc(5);

[cp,ct,dcp,dct] = cpvwb (cpct,R,Vinf,W0,0);
Pa0 = cp*0.5*dens*Area*(Vinf^3);

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
     0.05*(2*pi)                        * time; ...
     0.24*(2*pi)                        * time; ...
     0.04; ...
     0.40];

Ppcc = uc(1);
v = 0;
asum = 31;

dt = 0.05                               / time;
Nt = 100;
t = 0;
ts = dt*[0:Nt-1];

xs = [xt;xo;xc];
ys = [yt;Phat];

[dxdt,Phat,aac,bbc,ccc,ddc] = plantControl (xc,uc,p,atables,cpct,1);
[ACT,BCT,dt] = SSDiscreteTime (2,aac,bbc,dt);  % Discrete time.

Nxc = size(aac,1);
Nuc = size(bbc,2);
Nyc = size(ccc,1);

Nxs = Nxt + Nxo + Nxc;
Nys = Nyt+Nyo+Nyc;

% y vector:        u vector:
% (controls)       muV      1
% Phat     1       Dam      2
% Pe       2       a1sum    3
% beta     3       a2sum    4
% yaw      4       PhPCC    5
% (sensors)        PmPCC    6
% W        5       V        7  % Rotor-avg V at turbine.
% vnacx,y 6-7
% V,thV   8-9
% (estimates)
% Vst     10
% FTst    11
% Past    12
% (grid)
% ist   13,14
Ny = 14;
Nu = 7;
aa = zeros(Nxs,Nxs);
bbu = zeros(Nxs,Nu);
bby = zeros(Nxs,Ny);
cc = zeros(Ny,Nxs);
ddu = zeros(Ny,Nu);
ddy = zeros(Ny,Ny);


% Turbine and its connections to interface variables.
aa(1:Nxt,1:Nxt) = AT;
bbu(1:Nxt,7)  = BT(:,2);                         % External wind.
bby(1:Nxt,1)  = BT(:,1);                         % Power command from control.
cc(2:9,1:Nxt) = Cq([4, 2, 3, 1, 5, 6, 7, 8],:);  % Controls and sensors.
cc(13:14,1:Nxt) = Cq(60:61,:);                   % Grid currents.
ddu(2:9,7)    = Dq([4, 2, 3, 1, 5, 6, 7, 8],2); 
ddu(13:14,7)  = Dq(60:61,2);
ddy(2:9,1)    = Dq([4, 2, 3, 1, 5, 6, 7, 8],1);
ddy(13:14,1)  = Dq(60:61,1);

% Observer.
aa(Nxt+[1:Nxo],Nxt+[1:Nxo]) = AG;
% (no external inputs.)
bby(Nxt+[1:Nxo],2:4) = BG;  % Control inputs.
bby(Nxt+[1:Nxo],5:9) = GG;  % Sensor inputs.
cc(10:12,Nxt+[1:Nxo]) = CG([7, 10, 11],:);
ddy(10:12,2:4) = DG([7, 10, 11],:);

% Controller.
aa(Nxt+Nxo+[1:Nxc],Nxt+Nxo+[1:Nxc]) = ACT;
bbu(Nxt+Nxo+[1:Nxc],1:6) = BCT(:,[3, 7, 8, 9, 1, 2]);  % muV, Dam, a1sum, a2sum, PhPCC, PmPCC
bby(Nxt+Nxo+[1:Nxc],11:12) = BCT(:,[5, 4]);         % FTst, Past.
bby(Nxt+Nxo+[1:Nxc],5) = BCT(:,6);                  % Omega.
cc(1,Nxt+Nxo+[1:Nxc]) = ccc;
ddu(1,1:6) = ddc(:,[3, 7, 8, 9, 1, 2]);
ddy(1,11:12) = ddc(:,[5, 4]);
ddy(1,5) = ddc(:,6);

[ACL,BCL,CCL,DCL] = modularToUnifiedStateSpace(aa,bbu,bby,cc,ddu,ddy,ones(Ny,1));

eval(["save('-binary','ACL" inpnm ".bin','ACL');"]);
eval(["save('-binary','BCL" inpnm ".bin','BCL');"]);
eval(["save('-binary','CCL" inpnm ".bin','CCL');"]);
eval(["save('-binary','DCL" inpnm ".bin','DCL');"]);
