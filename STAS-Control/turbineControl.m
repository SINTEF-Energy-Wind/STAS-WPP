function [dxdt,yout,AA,BB,CC,DD,blydof,bludof] =                  ...
               turbineControl (linFlag,x,u,p,cpct,                ...
                               KeTab,WVTab,WPTab,bminTab,KTables, ...
                               KFTab,KSTab,KSqTab,KpiTab,KiiTab,RSCFlag)
%
% Build and link the wind-turbine-level controls.  This outputs pitch
% commands to the actuators, current commands to the generator
% converter control (part of the electric system), and the yaw angle
% command to the yaw actuator (not yet implemented).
%
%   States:              y vector:             u vector:
%
% --------------- Top-level inputs, outputs --------------
%                        bhat     1:3          Pc     1   (grid)
%                        ihgd,q   4:5          azi    2   (turbine)
%                        yhat      6           W      3   (turbine)
%                                              bet   4:6  (turbine)
%                                              Pe     7   (turbine)
%                                              vT    8:9  (turbine)
%                                              Mbl  10:12 (turbine)(*)
%                                              wang  13   (systems: wind angle for yaw.)
%
% ------------------ Windspeed observer ------------------
%   W*        1                                W      1  (u)
%   V*        2                                bet    2  (u, avg.)
%                                              Pem    3  (gen P)
%
% ------------------ Rotor speed control -----------------
%   Wm        1          bhRSC     1           W      1  (u)
%   betm      2          PehRSC    2           bet    2  (u, avg.)
%   PsiWb     3                                V*     3  (obs)
%   Pf        4                                Pc     4  (u)
%   PehRSC    5
%   ntw      6,7
%   ntb      8,9
%   blp      10
%
% -------------- Virtual induction generator -------------
%   wgm       1          Tgi(**)   1           wg     1  (u)
%   iwg       2
%
% ----------------- Active power control -----------------
%   Pem       1          igq       1           Phat   1  (multiple)
%   PsiP      2                                Pe     2  (u)
%   nd       3,4
%
% ------------------- Tower FA damping -------------------
%   vmF       1          bF        1           vF     1  (u)
%   ivmF      2                                betm   2  (RSC)
%   vdF       3
%
% ---------------- Tower SS damping (gen) ----------------
%   vmS       1          TgS       1           vS     1  (u)
%   ivmS      2                                V*     2  (obs)
%   vdS       3
%
% ------------------------ IBPSS -------------------------
%   vmS       1          betq     1:3          vS     1  (u)
%   ivmS      2                                azi    2  (u)
%   vdS       3                                betm   3  (RSC)
%
% ------------------------ IBPLF -------------------------
%   Mpsim    1:2         beti     1:3          M     1:3 (u)
%   PsiM     3:4                               azi    4  (u)
%                                              betm   5  (RSC)
%
% ------------------------- Yaw --------------------------
%
%
% (*) Root bending strain or nodal rotation may also be used, provided
% that the gains are scaled accordingly.
% (**) The gain is defined such that this is a current output that can
% be fed straight into the ihg,q output.
%
% Version:        Changes:
% --------        -------------
% 02.02.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 02.02.2019      Derivatives verified with complex step.
%
% Inputs:
% -------
% x               : 1: 
% u               : 1: 
% p               : 1: np    (-)       Number of generator poles.
%                   -------------------------
%                   1: Ro    (m)       Rotor outer radius.
%                   2: J     (kg m^2)  Rotor inertia.
%                   3: rho   (kg/m^3)  Density.
%                   4: KW    (-)       Gain on speed error.
%                   5: KV    (m)       Gain on speed error for windspeed.
%                   -------------------------
%                   1: aW    (rad/s)   LP filter on speed.
%                   2: ab    (rad/s)   LP filter on pitch.
%                   3: a1    (rad/s)   First filter on power.
%                   4: a2    (rad/s)   Second filter on power.
%                   5: Pr    (W)       Rated power.
%                   7: Wr    (rad/s)   Rated speed.
%                   6: fcd   (rad/s)   Rated speed fraction at transition to pitch control.
%                   8: anw   (rad/s)   Notch filter on speed.
%                   9: z1nw  (-)
%                  10: z2nw  (-)
%                  11: anb   (rad/s)   Notch filter on pitch.
%                  12: z1nb  (-)
%                  13: z2nb  (-)
%                  14: ablp  (rad/s)   Low-pass on pitch output.
%                   -------------------------
%                   1: ag    (rad/s)   BP filter frequency.
%                   2: zetag (-)       BP filter damping.
%                   3: Kd    (Nms/rad) Induction gen. stiffness.
%                   -------------------------
%                   1: Kp    (A/W)     Prop. gain on gen. power.
%                   2: Ki    (A/Ws)    Int. gain on gen. power.
%                   3: ap    (rad/s)   LP filter on power.
%                   4: anp   (rad/s)   Notch filter on power.
%                   5: z1np  (-)
%                   6: z2np  (-)
%                   -------------------------
%                   1: aF    (rad/s)   Frequency in rad/s.
%                   2: zetaF (-)       Damping ratio.
%                   3: azF   (rad/s)   Numerator frequency of phase shift.
%                   4: apF   (rad/s)   Denominator frequency of phase shift.
%                   -------------------------
%                   1: aS    (rad/s)   Frequency in rad/s (gen torque).
%                   2: zetaS (-)       Damping ratio.
%                   3: azS   (rad/s)   Numerator frequency of phase shift.
%                   4: apS   (rad/s)   Denominator frequency of phase shift.
%                   -------------------------
%                   1: aq    (rad/s)   Frequency in rad/s (IBP).
%                   2: zetaq (-)       Damping ratio.
%                   3: azq   (rad/s)   Numerator frequency of phase shift.
%                   4: apq   (rad/s)   Denominator frequency of phase shift.
%                   -------------------------
%                   1: aLP   (rad/s)   LP filter frequency.
%
% Outputs:
% --------
% 

thrd = 1/3;

% bl(y,u)dof: Blade 1, 2, and 3 DOFs for MBC transform.
blydof = [1:3];
bludof = [[4:6]; [10:12]];

Nu  = 13;   % Global input.
Ny  = 6;    % Global output.

Nxo = 2;
Nxr = 10;
Nxv = 2;
Nxg = 4;
Nxf = 3;
Nxs = 3;
Nxq = 3;
Nxi = 4;

Nuo = 3;  % Local u's.
Nur = 4;
Nuv = 1;
Nug = 2;
Nuf = 2;
Nus = 2;
Nuq = 3;
Nui = 5;

Nyo = 0;  % Local y's.
Nyr = 2;
Nyv = 1;
Nyg = 1;
Nyf = 1;
Nys = 1;
Nyq = 3;
Nyi = 3;

Npo = 5;
Npr = 14;
Npv = 3;
Npg = 6;
Npf = 4;
Nps = 4;
Npq = 4;
Npi = 1;

ixo = 0;
ixr = ixo + Nxo;
ixv = ixr + Nxr;
ixg = ixv + Nxv;
ixf = ixg + Nxg;
ixs = ixf + Nxf;
ixq = ixs + Nxs;
ixi = ixq + Nxq;

iuo = 0;
iur = iuo + Nuo;
iuv = iur + Nur;
iug = iuv + Nuv;
iuf = iug + Nug;
ius = iuf + Nuf;
iuq = ius + Nus;
iui = iuq + Nuq;

iyo = 0;
iyr = iyo + Nyo;
iyv = iyr + Nyr;
iyg = iyv + Nyv;
iyf = iyg + Nyg;
iys = iyf + Nyf;
iyq = iys + Nys;
iyi = iyq + Nyq;

np  = p(1);
ipo = 1;
ipr = ipo + Npo;
ipv = ipr + Npr;
ipg = ipv + Npv;
ipf = ipg + Npg;
ips = ipf + Npf;
ipq = ips + Nps;
ipi = ipq + Npq;

id = ['o';'r';'v';'g';'f';'s';'q';'i'];
Nid = size(id,1);
Nx = 0;
nnu = 0;
nny1 = 0;
for jj = 1:Nid
   eval (['Nx = Nx + Nx' id(jj) ';']);
   eval (['nnu = nnu + Nu' id(jj) ';']);
   eval (['nny1 = nny1 + Ny' id(jj) ';']);
end

% Make an overall y vector of the global outputs plus all the local
% inputs/outputs.
iyy = 0;
iyu = iyy + Ny;
iy1 = iyu + nnu;
nny = Ny + nnu + nny1; 

Pc   = u(1);
azi  = u(2);
W    = u(3);
bet1 = u(4);
bet2 = u(5);
bet3 = u(6);
Pe   = u(7);
vF   = u(8);
vS   = u(9);
Mb1  = u(10);
Mb2  = u(11);
Mb3  = u(12);
wang = u(13);
bavg = (1/3)*(bet1+bet2+bet3);

dxdt = zeros(Nx,1);
yout = zeros(Ny,1);

% Windspeed observer.
opar = p(ipo+[1:Npo]);
xo = x(ixo+[1:Nxo]);
uo = [W; bavg; x(ixg+1)];
[dxodt,aao,bbo] = shaftFunc (xo,uo,opar,cpct,KeTab);

% Rotor speed control.
rpar = p(ipr+[1:Npr]);
xr = x(ixr+[1:Nxr]);
ur = [W; bavg; x(ixo+2); Pc];
[dxrdt,yrout,aar,bbr,ccr,ddr] = RSC (xr,ur,rpar, ...
                                     WVTab,WPTab,bminTab,KTables,RSCFlag);

% Virtual induction generator.
vpar = p(ipv+[1:Npv]);
xv = x(ixv+[1:Nxv]);
uv = 0.5*np*W;
[dxvdt,Tgi,aav,bbv,ccv] = virtualInduction (xv,uv,vpar);

% Tower FA damping.
fpar = p(ipf+[1:Npf]);
xf = x(ixf+[1:Nxf]);
uf = [vF; x(ixr+2)];
[dxfdt,betF,aaf,bbf,ccf,ddf] = FADamping (xf,uf,fpar,KFTab);

% Tower SS damping (gen).
spar = p(ips+[1:Nps]);
xs = x(ixs+[1:Nxs]);
us = [vS; x(ixo+2)];
[dxsdt,TgS,aas,bbs,ccs,dds] = SSDamping (xs,us,spar,KSTab);

% Active power control.
gpar = p(ipg+[1:Npg]);
xg = x(ixg+[1:Nxg]);
ug = [x(ixr+5)+TgS*x(ixo+1); Pe];
[dxgdt,ygout,aag,bbg,ccg,ddg] = genPcontrol (xg,ug,gpar);

% Tower SS damping (IBP).
qpar = p(ipq+[1:Npq]);
xq = x(ixq+[1:Nxq]);
uq = [vS; azi; x(ixr+2)];
[dxqdt,betq,aaq,bbq,ccq,ddq] = IBPSS (xq,uq,qpar,KSqTab);

% Low-frequency individual blade pitch.
ipar = p(ipi+[1:Npi]);
xi = x(ixi+[1:Nxi]);
ui = [Mb1; Mb2; Mb3; azi; x(ixr+2)];
[dxidt,beti,aai,bbi,cci,ddi] = IBPLF (xi,ui,ipar,KpiTab,KiiTab);

% Yaw. (Not implemented.)


dxdt(ixo+[1:Nxo]) = dxodt;
dxdt(ixr+[1:Nxr]) = dxrdt;
dxdt(ixv+[1:Nxv]) = dxvdt;
dxdt(ixg+[1:Nxg]) = dxgdt;
dxdt(ixf+[1:Nxf]) = dxfdt;
dxdt(ixs+[1:Nxs]) = dxsdt;
dxdt(ixq+[1:Nxq]) = dxqdt;
dxdt(ixi+[1:Nxi]) = dxidt;

yout = [yrout(1)+betF+betq+beti;0;ygout-Tgi;0];

if (linFlag == 1)

   A  = spalloc (Nx,Nx,0.1*Nx^2);
   Bu = spalloc (Nx,Nu,0.05*Nx*Nu);
   By = spalloc (Nx,nny,0.05*Nx*nny);
   C  = spalloc (nny,Nx,0.05*nny*Nx);
   Du = spalloc (nny,Nu,0.06*nny*Nu);
   Dy = spalloc (nny,nny,0.01*nny*nny);

   % ---------- Windspeed observer -----------
   ir = ixo + [1:Nxo];
   ic = ixo + [1:Nxo];
   A(ir,ic) = aao;
   ic = iyu + iuo + [1:Nuo];
   By(ir,ic) = bbo;

   % ----------------- RSC -------------------
   ir = ixr + [1:Nxr];
   ic = ixr + [1:Nxr];
   A(ir,ic) = aar;
   ic = iyu + iur + [1:Nur];
   By(ir,ic) = bbr;
   ir = iy1 + iyr + [1:Nyr];
   ic = ixr + [1:Nxr];
   C(ir,ic) = ccr;
   ic = iyu + iur + [1:Nur];
   Dy(ir,ic) = ddr;

   % -------- Virtual induction gen. ---------
   ir = ixv + [1:Nxv];
   ic = ixv + [1:Nxv];
   A(ir,ic) = aav;
   ic = iyu + iuv + [1:Nuv];
   By(ir,ic) = bbv;
   ir = iy1 + iyv + [1:Nyv];
   ic = ixv + [1:Nxv];
   C(ir,ic) = ccv;

   % ---------- Active power control ---------
   ir = ixg + [1:Nxg];
   ic = ixg + [1:Nxg];
   A(ir,ic) = aag;
   ic = iyu + iug + [1:Nug];
   By(ir,ic) = bbg;
   ir = iy1 + iyg + [1:Nyg];
   ic = ixg + [1:Nxg];
   C(ir,ic) = ccg;
   ic = iyu + iug + [1:Nug];
   Dy(ir,ic) = ddg;

   % ------------ Tower FA damping -----------
   ir = ixf + [1:Nxf];
   ic = ixf + [1:Nxf];
   A(ir,ic) = aaf;
   ic = iyu + iuf + [1:Nuf];
   By(ir,ic) = bbf;
   ir = iy1 + iyf + [1:Nyf];
   ic = ixf + [1:Nxf];
   C(ir,ic) = ccf;
   ic = iyu + iuf + [1:Nuf];
   Dy(ir,ic) = ddf;

   % -------- Tower SS damping (gen) ---------
   ir = ixs + [1:Nxs];
   ic = ixs + [1:Nxs];
   A(ir,ic) = aas;
   ic = iyu + ius + [1:Nus];
   By(ir,ic) = bbs;
   ir = iy1 + iys + [1:Nys];
   ic = ixs + [1:Nxs];
   C(ir,ic) = ccs;
   ic = iyu + ius + [1:Nus];
   Dy(ir,ic) = dds;

   % ---------------- IBPSS ------------------
   ir = ixq + [1:Nxq];
   ic = ixq + [1:Nxq];
   A(ir,ic) = aaq;
   ic = iyu + iuq + [1:Nuq];
   By(ir,ic) = bbq;
   ir = iy1 + iyq + [1:Nyq];
   ic = ixq + [1:Nxq];
   C(ir,ic) = ccq;
   ic = iyu + iuq + [1:Nuq];
   Dy(ir,ic) = ddq;

   % ---------------- IBPLF ------------------
   ir = ixi + [1:Nxi];
   ic = ixi + [1:Nxi];
   A(ir,ic) = aai;
   ic = iyu + iui + [1:Nui];
   By(ir,ic) = bbi;
   ir = iy1 + iyi + [1:Nyi];
   ic = ixi + [1:Nxi];
   C(ir,ic) = cci;
   ic = iyu + iui + [1:Nui];
   Dy(ir,ic) = ddi;

   % Link local u's to global u's.
   ir = iyu + [iuo+[1 2 2 2] iur+[1 2 2 2 4] iuv+1 iug+2 ...
               iuf+1 ius+1 iuq+[1 2] iui+[1 2 3 4]];
   ic = [3 4 5 6 3 4 5 6 1 3 7 8 9 9 2 10 11 12 2];
   ss = [1 thrd thrd thrd 1 thrd thrd thrd 1 0.5*np ...
         1 1 1 1 1 1 1 1 1];
   Du = Du + sparse(ir,ic,ss,nny,Nu);

   % Link local u's to x's.
   ir = iyu + [iuo+3 iur+3 iug+1 iuf+2 ius+2 iuq+3 iui+5];
   ic =       [ixg+1 ixo+2 ixr+5 ixr+2 ixo+2 ixr+2 ixr+2];
   ss = [1 1 1 1 1 1 1];
   C = C + sparse(ir,ic,ss,nny,Nx);

   % Link gen u's to local y's and x's.
   ir = iyu + iug + 1;
   ic = iy1 + iys + 1;
   Dy(ir,ic) = Dy(ir,ic) + x(ixo+1);
   ir = iyu + iug + 1;
   ic = ixo + 1;
   C(ir,ic) = C(ir,ic) + TgS;

   % Link global y's to local y's.
   ir = iyy + [1 1 1 1 2 2 2 2 3 3 3 3 5 5];
   ic = iy1 + [iyr+1 iyf+1 iyq+1 iyi+1 ...
               iyr+1 iyf+1 iyq+2 iyi+2 ...
               iyr+1 iyf+1 iyq+3 iyi+3 ...
               iyg+1 iyv+1];
   ss = [1 1 1 1 1 1 1 1 1 1 1 1 1 -1];
   Dy = Dy + sparse(ir,ic,ss,nny,nny);
   
   Dnorm = ones(nny,1);
   [AA,BB,CC,DD] = modularToUnifiedStateSpace (A,Bu,By,C,Du,Dy,Dnorm);

else

   AA = sparse (Nx,Nx);
   BB = sparse (Nx,Nu);
   CC = sparse (nny,Nx);
   DD = sparse (nny,Nu);

end


