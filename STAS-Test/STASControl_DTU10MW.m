function c = STASControl_DTU10MW ()

load 'LTMnorms.txt';
[length,time,mass,current,voltage,        ...
 velocity,force,power,stress,ndens,nvisc, ...
 stiffness,damping,resistance,inductance, ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

load 'cpct_DTU10MW.txt';
c.cpct = cpct_DTU10MW;


c.RSCFlag = 0;


% These are the parameters used in the controller.
% -----------------------
c.np    = 198;
% -----------------------
c.Ro    = 89.15                   / length;
c.J     = 1.6e8                   / (mass/(length^2));
c.rho   = 1.225                   / (mass/(length^3));
c.KW    = 2;
c.KV    = 15                     / length;
% -----------------------
c.aW    = 0.4*(2*pi)              * time;
c.ab    = 0.5*(2*pi)              * time;
c.a1    = 0.52*(2*pi)              * time;
c.a2    = 0.48*(2*pi)              * time;
c.Pr    = 1e7                     / power;
c.Wr    = 1.0053                  * time;
c.fcd   = 0.95;
c.anw   = 0.238*(2*pi)             * time;  % CHECK
c.z1nw  = 0.01;
c.z2nw  = 0.10;
c.anb   = 0.239*(2*pi)             * time;  % CHECK
c.z1nb  = 0.04;
c.z2nb  = 0.40;
c.ablp  = 2*(2*pi)                * time;  % CHECK
% -----------------------
c.ag    = 1.9*(2*pi)              * time;   % CHECK
c.zetag = 0.1;
c.Kd    = 2e7                     / (force*length*time);   % CHECK, units for current output.
% -----------------------
c.Kp    = -0.00015                / (current/power);
c.Ki    = -0.0025                 / (current/(power*time));
c.ap    = 1.5*(2*pi)              * time;   % CHECK
c.anp   = c.ag;
c.z1np  = 0.01;
c.z2np  = 0.10;
% -----------------------
c.aF    = 0.238*(2*pi)            * time;   % CHECK
c.zetaF = 0.1;
c.azF   = 0.2*(2*pi)              * time;   % CHECK
c.apF   = 0.4*(2*pi)              * time;   % CHECK
% -----------------------
c.aS    = 0.239*(2*pi)            * time;   % CHECK
c.zetaS = 0.1;
c.azS   = 0.22*(2*pi)             * time;   % CHECK
c.apS   = 0.28*(2*pi)             * time;   % CHECK
% -----------------------
c.aq    = 0.237*(2*pi)            * time;   % CHECK
c.zetaq = 0.1;                              % CHECK
c.azq   = 0.2*(2*pi)              * time;   % CHECK
c.apq   = 0.2*(2*pi)              * time;   % CHECK
% -----------------------
c.aLP   = 0.08*(2*pi)             * time;   % CHECK

% -----------------------------------------------------------------------
% Other parameters, used below, that might also be useful to output.
c.Wmin  = 0.628                   *time;
c.fab   = 1.1;
c.fbc   = 0.9;
c.dPdWr = 1e5                     /(power*time);
c.TSR   = 7.5;
c.Vhdc  = 6500                    / voltage;

c.cpar = [c.np;c.Ro;c.J;c.rho;c.KW;c.KV;c.aW;c.ab;c.a1;c.a2;c.Pr;c.fcd;c.Wr; ...
          c.anw;c.z1nw;c.z2nw;c.anb;c.z1nb;c.z2nb;c.ablp;c.ag;c.zetag;c.Kd;  ...
          c.Kp;c.Ki;c.ap;c.anp;c.z1np;c.z2np;c.aF;c.zetaF;c.azF;c.apF;       ...
          c.aS;c.zetaS;c.azS;c.apS;c.aq;c.zetaq;c.azq;c.apq;c.aLP];

Area  = pi*(c.Ro^2)               /(length^2);

% -----------------------------------------------------------------------
% KeTab: Table of electrical efficiency as a function of power.  Used in
% the windspeed observer.
c.KeTab = [([[0:0.2:1] [2:20]].')*1e6/power                       ...
           [0.600 0.720 0.800 0.870 0.890 0.909 0.930 0.935 0.940 ...
            0.940 0.940 0.938 0.938 0.936 0.934 0.933 0.930 0.925 ...
            0.921 0.917 0.914 0.910 0.906 0.901 0.896].'];

% -----------------------------------------------------------------------
% WVTab: table of rotor speed as a function of windspeed.  Used in RSC
% during curtailed operation.
Nb = 10;
Wbs = linspace(c.fab*c.Wmin,c.fbc*c.Wr,Nb);
Vbs = c.Ro*Wbs/c.TSR;

c.WVTab = zeros(30,2);
c.WVTab(:,1) = [0 [3:30] 50].'     / velocity;
c.WVTab(:,2) = max(min(1.02*c.TSR*c.WVTab(:,1)/c.Ro,c.Wr),c.fab*c.Wmin);

% -----------------------------------------------------------------------
% WPTab: The primary schedule of electrical power as a function of 
% rotor speed.
etag = 0.925;   % TEMPORARY.
c.WPTab = [[0 c.Wmin Wbs c.fcd*c.Wr c.Wr 2*c.Wr].'                      ...
           [-c.dPdWr*c.Wmin, 0, etag*0.45*0.5*c.rho*Area*(Vbs.^3), ...
            c.Pr+c.dPdWr*(c.fcd*c.Wr-c.Wr), c.Pr, c.Pr+c.dPdWr*c.Wr].'];

% -----------------------------------------------------------------------
% bminTab: The minimum pitch angle as a function of windspeed.
c.bminTab = [[0 [3:30] 50].'/velocity ...
              (pi/180)*[3 3 2.751 1.966 0.896 0 0 0 0.5 2.0 3*ones(1,19) 3].'];

% -----------------------------------------------------------------------
% KpbTab, KibTab, KdbTab: Blade pitch proportional, integral, and derivative
% gain schedules as a function of pitch angle.
bref = 23;
Nk = bref + 3;
KpbTab = zeros(Nk,2);
KpbTab(:,1) = [(pi/180)*[-90 [0:1:bref] 90]].';
KpbTab(:,2) = 1 + KpbTab(:,1).*(-2.541 + KpbTab(:,1).*(-7.814 + KpbTab(:,1).* ...
              (46.281 + KpbTab(:,1).*(-59.871)))); 
KpbTab(1,2)  = KpbTab(2,2) - 0.01*(KpbTab(1,1)  - KpbTab(2,1));
KpbTab(Nk,2) = KpbTab(Nk-1,2) - 0.01*(KpbTab(Nk,1) - KpbTab(Nk-1,1));
KpbTab(:,2) = KpbTab(:,2)           / time;

KibTab = zeros(Nk,2);
KibTab(:,1) = [(pi/180)*[-90 [0:1:bref] 90]].';
KibTab(:,2) = 0.351 + KibTab(:,1).*(-2.405 + KibTab(:,1).*(13.128 + KibTab(:,1).* ...
              (-31.926 + KibTab(:,1).*(27.689))));
KibTab(1,2)  = KibTab(2,2) - 0.01*(KibTab(1,1)  - KibTab(2,1));
KibTab(Nk,2) = KibTab(Nk-1,2) - 0.01*(KibTab(Nk,1) - KibTab(Nk-1,1));
KibTab(:,2) = KibTab(:,2)           / (time^2);

KdbTab = zeros(Nk,2);
KdbTab(:,1) = [(pi/180)*[-90 [0:1:bref] 90]].';
KdbTab(:,2) = KpbTab(:,2);
KdbTab(1,2)  = KdbTab(2,2) - 0.01*(KdbTab(1,1)  - KdbTab(2,1));
KdbTab(Nk,2) = KdbTab(Nk-1,2) - 0.01*(KdbTab(Nk,1) - KdbTab(Nk-1,1));

% For quick scaling of the pitch gains...
%-------------------------
KpbTab(:,2) = 1*KpbTab(:,2);
KibTab(:,2) = 1*KibTab(:,2);
KdbTab(:,2) = 0*KdbTab(:,2);   % Deactivate derivative gain.
%-------------------------

% -----------------------------------------------------------------------
% KFTab: Gains on tower fore-aft damping control, as a function of blade
% pitch angle.
c.KFTab = zeros(Nk,2);
c.KFTab(:,1) = [(pi/180)*[-90 [0:1:bref] 90]].';
c.KFTab(:,2) = 0.005                           * length;    % NEEDS SCHEDULING.

% -----------------------------------------------------------------------
% KSTab: Gains on tower side-to-side damping control using the generator,
% as a function of windspeed.
c.KSTab = zeros(30,2);
c.KSTab(:,1) = [0 [3:30] 50].'                 / velocity;
c.KSTab(:,2) = 5e6                             / force;     % NEEDS SCHEDULING.

% -----------------------------------------------------------------------
% KSqTab: Gains on tower side-to-side damping control using individual
% blade pitch, as a function of pitch.
c.KSqTab = zeros(Nk,2);
c.KSqTab(:,1) = [(pi/180)*[-90 [0:1:bref] 90]].';
c.KSqTab(:,2) = 0.005                          * velocity;  % CHECK.

% -----------------------------------------------------------------------
% KpiTab, KiiTab: Gains on individual pitch control based on blade root
% moments, strain, or nodal rotation.  Scheduled as a function of pitch.
c.KpiTab = zeros(Nk,2);
c.KpiTab(:,1) = [(pi/180)*[-90 [0:1:bref] 90]].';
c.KpiTab(:,2) = 40;                                         % CHECK.  Based on nodal rotation.

c.KiiTab = zeros(Nk,2);
c.KiiTab(:,1) = [(pi/180)*[-90 [0:1:bref] 90]].';
c.KiiTab(:,2) = 30;                                         % CHECK.  Based on nodal rotation.

% -----------------------------------------------------------------------
% Plain tables.
c.KpbTab0 = KpbTab;
c.KibTab0 = KibTab;
c.KdbTab0 = KdbTab;
c.KeTab0 = c.KeTab;
c.WVTab0 = c.WVTab;
c.WPTab0 = c.WPTab;
c.bminTab0 = c.bminTab;
c.KFTab0 = c.KFTab;
c.KSTab0 = c.KSTab;
c.KSqTab0 = c.KSqTab;
c.KpiTab0 = c.KpiTab;
c.KiiTab0 = c.KiiTab;

% -----------------------------------------------------------------------
% pchip spline tables.
c.KeTab   = interp1 (c.KeTab(:,1),  c.KeTab(:,2),  'pchip','pp');
c.WVTab   = interp1 (c.WVTab(:,1),  c.WVTab(:,2),  'pchip','pp');
c.WPTab   = interp1 (c.WPTab(:,1),  c.WPTab(:,2),  'pchip','pp');
c.bminTab = interp1 (c.bminTab(:,1),c.bminTab(:,2),'pchip','pp');

KpbTab  = interp1 (KpbTab(:,1), KpbTab(:,2), 'pchip','pp');
KibTab  = interp1 (KibTab(:,1), KibTab(:,2), 'pchip','pp');
KdbTab  = interp1 (KdbTab(:,1), KdbTab(:,2), 'pchip','pp');
c.KTables = [KpbTab KibTab KdbTab];

c.KFTab   = interp1 (c.KFTab(:,1),  c.KFTab(:,2),  'pchip','pp');
c.KSTab   = interp1 (c.KSTab(:,1),  c.KSTab(:,2),  'pchip','pp');
c.KSqTab  = interp1 (c.KSqTab(:,1), c.KSqTab(:,2), 'pchip','pp');
c.KpiTab  = interp1 (c.KpiTab(:,1), c.KpiTab(:,2), 'pchip','pp');
c.KiiTab  = interp1 (c.KiiTab(:,1), c.KiiTab(:,2), 'pchip','pp');

