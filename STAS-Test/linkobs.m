
clear;

[length,time,mass,current,voltage,              ...
 velocity,force,power,stress,density,viscosity, ...
 stiffness,damping,resistance,inductance,       ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

nm = 'DTU10MW';
inpnm = '_P060_';
outnm = '_P060_V100';

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

KI = diag([0.01, 0.1, 0.1, 0.1, 0.1]);
eval(["save('-binary','KI" outnm ".bin','KI');"]);

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

eval(["load 'AT" txt Vstr ".bin';"]);
eval(["load 'BT" txt Vstr ".bin';"]);
eval(["load 'Cq" txt Vstr ".bin';"]);
eval(["load 'Dq" txt Vstr ".bin';"]);
eval(["load 'AG" txt Vstr ".bin';"]);
eval(["load 'BG" txt Vstr ".bin';"]);
eval(["load 'CG" txt Vstr ".bin';"]);
eval(["load 'DG" txt Vstr ".bin';"]);
eval(["load 'GG" txt Vstr ".bin';"]);

% u(inputs) : 1: Pc, 2: V, 3: thV, 4: Fw, 5-6: vs.
% y(sensors): 1: W, 2: beta0, 3: yaw, 4: Pe, 5-6: vnac, 7: V, 8: thV
% y(other)  : 9: a, 10: FT, 11: Pa, 12-23: qfnd, 24-59: qbl, 60-61: isd,q.
iyc = [4, 2, 3].';
iys = [1, 5, 6, 7, 8].';

Nxt = size(AT,1);
Nxo = size(AG,1);
Ny  = size(Cq,1);
Nu  = size(BT,2);
Nys = size(GG,2);
Nyc = size(iyc,1);

% Check some step functions.
Nt = 2400;
yt = zeros(Nt,Ny*Nu);
yo = zeros(Nt,Ny*Nu);
xt = zeros(Nxt,Nu);
xo = zeros(Nxo,Nu);
xi = zeros(Nys,Nu);
u  = zeros(Nu,1);
ts = zeros(Nt,1);
for it = 2:Nt

   ts(it) = ts(it-1) + dt;

   if (it > 200)
      u = ones(Nu,1);
   end

%   if (it < 30/dt)
%      u = (it*dt/30)*ones(Nu,1);
%   end

%   u = sin(2*pi*0.01*ts(it))*ones(Nu,1);

   for iu = 1:Nu
      icn = Ny*(iu-1);
      u1 = zeros(Nu,1);
      u1(iu) = u(iu);
      xt(:,iu) = AT*xt(:,iu) + BT*u1;
      yt(it,icn+[1:Ny]) = (Cq*xt(:,iu) + Dq*u1).';

      % Pass to the turbine: sensor measurements and control signals.
      uc = yt(it,icn+iyc).';  % Controls: P, beta, yaw.

      xo(:,iu) = AG*xo(:,iu) + BG*uc + GG*(yt(it,icn+iys).' + xi(:,iu));
      yo(it,icn+[1:Ny]) = (CG*xo(:,iu) + DG*uc).';
      xi(:,iu) = xi(:,iu) + KI*(yt(it,icn+iys).' - yo(it,icn+iys).');
   end

end

iyplt = [9, 10, 11, 60, 15, 16, 34, 35];
posz = [0.100 0.810 0.355 0.150; ...
        0.100 0.620 0.355 0.150; ...
        0.100 0.430 0.355 0.150; ...
        0.100 0.240 0.355 0.150; ...
        0.55 0.810 0.355 0.150; ...
        0.55 0.620 0.355 0.150; ...
        0.55 0.430 0.355 0.150; ...
        0.55 0.240 0.355 0.150];
for iu = 1:2 % 6
   icn = Ny*(iu-1);
   figure(iu);
   clf;
   set (gcf,'papersize',[8. 11.],'paperorientation','portrait', ...
        'paperposition',[0.1 0.2 8. 11.]);
   hax = zeros(8,1);
   for iyout = 1:8
      if (iyout == 4) || (iyout == 8)
         xtck = [0:20:120];
      else
         xtck = [];
      end
      hax(iyout) = axes();
      set(hax(iyout),'position',posz(iyout,:),'xlim',[0 120], ...
          'fontsize',10,'fontname','Times-Roman','ticklength',[0 0], ...
          'xtick',xtck,'linewidth',1);
      axes(hax(iyout));
      hold on;
      box on;
      plot (ts,yt(:,icn+iyout),'color',[0.0, 0.0, 0.0],'linewidth',2.0);
      plot (ts,yo(:,icn+iyout),'color',[0.6, 0.6, 0.6],'linewidth',2.0);
      hold off;
   end
   eval(["saveas (gcf,'step" int2str(iu) ".pdf','pdf');"]);
   figure(6+iu);
   clf;
   set (gcf,'papersize',[8. 11.],'paperorientation','portrait', ...
        'paperposition',[0.1 0.2 8. 11.]);
   for iyout = 1:8
      hax(iyout) = axes();
      set(hax(iyout),'position',posz(iyout,:),'xlim',[0 120], ...
          'fontsize',10,'fontname','Times-Roman','ticklength',[0 0], ...
          'xtick',[0:20:120],'linewidth',1);
      axes(hax(iyout));
      hold on;
      box on;
      plot (ts,yt(:,icn+iyplt(iyout)),'color',[0.0, 0.0, 0.0],'linewidth',2.0);
      plot (ts,yo(:,icn+iyplt(iyout)),'color',[0.6, 0.6, 0.6],'linewidth',2.0);
      hold off;
   end
   eval(["saveas (gcf,'step" int2str(4+iu) ".pdf','pdf');"]);
end

