function [s,a] = STASTurbine_NREL5MW_land ()

% -----------------------------------------------
% NREL 5 MW Reference Wind Turbine, land-based tower,
% free foundation (no soil restraints).
% -----------------------------------------------


%==========================================================
% Normalization.
%==========================================================
length = 1
time   = 1
power  = 1e6
voltage = sqrt(power)
%---------------------------------------
velocity = length/time
mass   = power*(time^3)/(length^2)
force  = mass*length/(time^2)
stress = force/(length^2)
ndens  = mass/(length^3)
nvisc  = mass/(length*time)
stiffness = force/length
damping = force*time/length
current = power/voltage
resistance = voltage/current
inductance = voltage*time/current
capacitance = current*time/voltage
flux   = voltage*time

LTMnorms = [length;time;mass;current];
save('-ascii','LTMnorms.txt','LTMnorms');

%==========================================================
% General.
%==========================================================
Nef = 20;    % Check that Nef is sufficient for wave loads.
Nmud = 5;    % (Nmud is part of Nef, specifies number below seabed.)
Nwater = 10; % (Nwater  "      "   , number between seabed and waterline.)
Net = 10;
Nev = 2;
Ner = 2;
Nen = 4;
Ned = 4;
Neh = 2;
Neb = 16;

% Number of modes to retain.  Retain as few as possible in order
% that the numerical conditioning of the structural matrices
% will be as good as possible.  (This can also be accomplished
% by later partitioning of the system matrices, if, for instance,
% a study is focused on a group of high-frequency modes.)
%              |  Caution, if you use more than this number, the
%              V  code has not been validated.
Nfnd = 10; %   20;  % Foundation.
Ntow = 10; %   20;  % Tower.
Nnac = 2;  %   8;   % Nacelle.
Ndrv = 2;  %   18;  % Driveshaft.
Nbld = 16; %   20;  % Blades.

Nwel = 15; % Number of elements subject to wave loads.

% Set to positive integers if icp contains the control points
% for a spline representation of the aerodynamic states. Set
% to negative integers if icp contains the blade body modes
% for a modal-slave representation of the aerodynamic states.
% In the latter case, the selection of modes has to be verified
% as the turbine design is altered.
%icp = round(1 + (Neb-1)*[0.00 0.25 0.50 0.70 0.85 1.00]'); % Splines.
%icp = -[1 3 5 7 8 9]'; % Modes. Verify as design is altered.
icp = -1;
%icp = 0;

% Initialization.
s = createStructure (Nef,Net,Nev,Ner,Nen,Ned,Neh,Neb,Nwel);

s.zdamp = 0.008;      % Modal structural damping ratio.
                      % or ...
s.adamp = 0;          % Mass-proportional damping factor.
s.bdamp = 0; % 0.0002;     % Stiffness-proportional damping factor.

% Correct zero-lift aoa values for inboard DTU 10 MW profiles.
load('-binary','aoazs_NREL5MW.bin');
Nfoils = size(aoazs,1);
a = createAero (Nfoils,Neb);

s.foundation.Nmod   = Nfnd;
s.tower.Nmod        = Ntow;
s.nacelle.Nmod      = Nnac;
s.driveshaft.Nmod   = Ndrv;
for ib = 1:s.Nb
   s.blade(ib).Nmod = Nbld;
end

%==========================================================
% Blades.
%==========================================================

Lb = 61.5   /length; % Blade length from pitch bearing flange to tip.

% Section properties.  In the structural calculations the structural
% twist defines the section coordinate system, while in the
% aerodynamic calculations the aerodynamic twist defines the section
% coordinate system.
% 
% Input properties associated with each node in the blade.  There are
% Neb+1 nodes.
%
% x/Lb egy xis(deg) rhoA rhoJ EA EIyy EIzz GJ
props = [ ...
0.0000E+00  -0.00017  1.3308E+01  6.7894E+02  5.0080E+03  9.7295E+09  1.8110E+10  1.8114E+10  5.5644E+09; ...
5.2030E-02  -0.04345  1.3308E+01  7.4004E+02  4.1999E+03  9.8678E+09  1.5287E+10  1.9789E+10  4.6666E+09; ...
1.1707E-01  -0.08394  1.3308E+01  4.0064E+02  1.4133E+03  4.4940E+09  5.5284E+09  8.0632E+09  1.5704E+09; ...
1.8211E-01  -0.22235  1.3181E+01  4.1682E+02  6.0504E+02  4.0824E+09  3.9495E+09  7.2717E+09  6.7227E+08; ...
2.4715E-01  -0.25941  1.1072E+01  3.4948E+02  2.8022E+02  3.0116E+09  2.3887E+09  4.9485E+09  3.1135E+08; ...
3.2846E-01  -0.20382  9.6720E+00  3.3000E+02  2.0594E+02  2.3575E+09  1.8283E+09  4.2441E+09  2.2882E+08; ...
3.9350E-01  -0.19323  8.5340E+00  3.1382E+02  1.5694E+02  1.9441E+09  1.3619E+09  3.7508E+09  1.7438E+08; ...
4.9106E-01  -0.13252  6.7110E+00  2.6334E+02  7.3071E+01  1.1688E+09  6.8130E+08  2.7342E+09  8.1190E+07; ...
5.5610E-01  -0.14035  5.5460E+00  2.4167E+02  5.1705E+01  9.2295E+08  4.0890E+08  2.3340E+09  5.7450E+07; ...
6.5366E-01  -0.17418  3.8340E+00  1.7940E+02  2.4696E+01  5.3970E+08  1.7588E+08  1.3234E+09  2.7440E+07; ...
7.1870E-01  -0.26022  2.8900E+00  1.5441E+02  1.6686E+01  4.6001E+08  1.0726E+08  1.0202E+09  1.8540E+07; ...
7.8376E-01  -0.22795  2.1160E+00  1.2956E+02  1.3077E+01  3.2889E+08  7.6310E+07  7.0961E+08  1.4530E+07; ...
8.4878E-01  -0.21662  1.3420E+00  9.8776E+01  7.2540E+00  2.1160E+08  4.9480E+07  4.5487E+08  8.0600E+06; ...
8.9756E-01  -0.23124  7.6000E-01  8.3001E+01  5.4810E+00  1.6025E+08  3.4670E+07  3.5372E+08  6.0900E+06; ...
9.4636E-01  -0.09470  2.5300E-01  5.9340E+01  3.8160E+00  6.3230E+07  1.9630E+07  1.5881E+08  4.2400E+06; ...
9.7886E-01  -0.07096  1.0100E-01  4.5818E+01  1.9530E+00  2.9920E+07  7.5500E+06  8.5070E+07  2.1700E+06; ...
1.0000E+00  -0.05181  0.0000E+00  1.0319E+01  1.7100E-01  3.5300E+06  1.7000E+05  5.0100E+06  1.9000E+05];

rbnod   = Lb*props(:,1)     /length;
rbel    = 0.5*(rbnod(2:Neb+1) + rbnod(1:Neb));

EA      =  props(:,6)/force;
EIyy    =  props(:,7)/(force*length^2);
EIzz    =  props(:,8)/(force*length^2);
GJ      =  props(:,9)/(force*length^2);
egy     =  props(:,2)/length;        % Centroid in structural section CS.
egz     =  zeros(Neb,1);
xis     =  zeros(Neb,1);
rhosA   =  props(:,4)/(mass/length); % Mass per unit length.
rhosJ   =  props(:,5)/(mass*length); % Torsional inertia per unit length.
cgy     =  zeros(Neb,1);
cgz     =  zeros(Neb,1);

% Aerodynamic properties.
% xia(deg), chord, xpc(pitch axis), t/c
propsn = [ ...
1.3308E+01  3.5420E+00  5.0000E-01  1.00; ...
1.3308E+01  3.8540E+00  4.9000E-01  1.00; ...
1.3308E+01  4.1670E+00  4.6000E-01  1.00; ...
1.2885E+01  4.5790E+00  4.1000E-01  0.40; ...
1.1207E+01  4.6118E+00  3.7000E-01  0.35; ...
9.6706E+00  4.3688E+00  3.6000E-01  0.35; ...
8.5216E+00  4.1516E+00  3.5000E-01  0.30; ...
6.7118E+00  3.7827E+00  3.5000E-01  0.25; ...
5.5485E+00  3.5410E+00  3.5000E-01  0.21; ...
3.8639E+00  3.1810E+00  3.5000E-01  0.21; ...
2.8989E+00  2.9410E+00  3.5000E-01  0.18; ...
2.1157E+00  2.7009E+00  3.5000E-01  0.18; ...
1.3417E+00  2.4610E+00  3.5000E-01  0.18; ...
7.6682E-01  2.2687E+00  3.5000E-01  0.18; ...
2.9262E-01  1.8905E+00  3.5000E-01  0.18; ...
1.0600E-01  1.4190E+00  3.5000E-01  0.18; ...
1.0600E-01  1.4190E+00  3.5000E-01  0.18];

props = 0.5*(propsn(1:Neb,:) + propsn(2:Neb+1,:));

chord   = props(:,2)/length;
xia     = props(:,1)*pi/180;     % Aero twist, rad.  Positive is LE into wind.
xpc     = props(:,3);            % Origin of section CS aft of the LE,
                                 % fraction of chord.
tca     = props(:,4);            % t/c.

% Element coordinates, pitch CS.
xbaero  = [rbel zeros(Neb,1) zeros(Neb,1)].';

% Aero center aft of LE, fraction of chord.
acent   = [ones(2,1);0.25*ones(14,1)];

% Nodal coordinates, pitch CS.
xbnaero = [rbnod zeros(Neb+1,1) zeros(Neb+1,1)].'; 
Lel     = norm(xbnaero(:,2:Neb+1) - xbnaero(:,1:Neb),2,'cols').';

% Compute foil weights from t/c.
load('-ascii','tcfoils_NREL5MW.txt');
foilwt = getFoilWeights (tcfoils_NREL5MW,tca);

% Put values into the data structures.
for ib = 1:s.Nb

   s.blade(ib).Lb      = Lb;
   s.blade(ib).Lel(:)  = Lel;
   s.blade(ib).xis(:)  = xis;
   s.blade(ib).zeta    = s.zdamp;

   % Store the element twist values in the outboard node.  This is because
   % Node 1 will take on a different meaning as the reference node.
   idof = [10:6:6*(Neb+1)-2].';
   s.blade(ib).Pn_B(idof)    = -xia + xis;             % CHECK SIGN ON XIS.
   s.blade(ib).Pn_B(4)       = -xia(1) + xis(1);

%   % Put the cone angle in the reference node part of Pn_B. This
%   % way, I don't need to pass around the cone angle as a separate
%   % parameter.
%   s.blade(ib).Pn_B(5) = -phi;

   s.blade(ib).Pn_B(1:3) = zeros(3,1);  % Make sure this is exactly zero.
   for inod = 2:Neb+1
      idof = 6*(inod-1);
      s.blade(ib).Pn_B(idof+1:idof+3) = xbnaero(:,inod);
   end

   s.blade(ib).conn = [ones(1,s.blade(ib).Nel); ...
                      [1:s.blade(ib).Nel];      ...
                      [2:s.blade(ib).Nel+1]];

   for iel = 1:Neb

      ic6   = 6*(iel-1);
      ic12  = 12*(iel-1);
      ic144 = 144*(iel-1);
      s.blade(ib).rhos(:,ic6+[1:6]) = ...
             diag([rhosA(iel);rhosA(iel);rhosA(iel);rhosJ(iel);0;0]);

      % Account for the center-of-mass offset.
      val = rhosA(iel)*cgy(iel);
      s.blade(ib).rhos(1,ic6+6) = -val;
      s.blade(ib).rhos(3,ic6+4) =  val;
      s.blade(ib).rhos(4,ic6+3) =  val;
      s.blade(ib).rhos(6,ic6+1) = -val;
      s.blade(ib).rhos(4,ic6+4) =  s.blade(ib).rhos(4,ic6+4) ...
                                +  rhosA(iel)*(cgy(iel)^2);

      val = rhosA(iel)*cgz(iel);
      s.blade(ib).rhos(1,ic6+5) =  val;
      s.blade(ib).rhos(2,ic6+4) = -val;
      s.blade(ib).rhos(4,ic6+2) = -val;
      s.blade(ib).rhos(5,ic6+1) =  val;
      s.blade(ib).rhos(4,ic6+4) =  s.blade(ib).rhos(4,ic6+4) ...
                                +  rhosA(iel)*(cgz(iel)^2);

      s.blade(ib).EEs(:,ic6+[1:6]) = ...
             diag([EA(iel);0;0;GJ(iel);EIyy(iel);EIzz(iel)]);

      % Account for the centroid offset.
      s.blade(ib).EEs(1,ic6+5)  =  EA(iel)*egz(iel);
      s.blade(ib).EEs(1,ic6+6)  = -EA(iel)*egy(iel);
      s.blade(ib).EEs(5,ic6+1)  =  EA(iel)*egz(iel);
      s.blade(ib).EEs(6,ic6+1)  = -EA(iel)*egy(iel);

      s.blade(ib).ke_s(:,ic12+[1:12]) = ...
             buildkes (s.blade(ib).EEs(:,ic6+[1:6]),s.blade(ib).Lel(iel));
      s.blade(ib).me_s(:,ic12+[1:12]) = ...
             buildmes (s.blade(ib).rhos(:,ic6+[1:6]),s.blade(ib).Lel(iel));

   end
end

a.chord(:)    = [chord;chord;chord];
a.xia(:)      = [xia;xia;xia];
a.xpc(:)      = [xpc;xpc;xpc];
a.acent(:)    = [acent;acent;acent];
a.foilwt(:,:) = [foilwt foilwt foilwt];
a.Lel(:)      = [Lel;Lel;Lel];
a.icp         = icp;

load('-binary','aoas_NREL5MW.bin');
load('-binary','kfoils_NREL5MW.bin');
a.aoas        = aoas;
a.kfoils      = kfoils;

a.Pa_B(4:6:6*a.Nb*Neb-2,1) = -a.xia;

for ib = 1:a.Nb
   for iel = 1:Neb
      idof = 6*Neb*(ib-1) + 6*(iel-1);
      ic   = 3*Neb*(ib-1) + 3*(iel-1);

      a.Pa_B(idof+1:idof+3,1) = xbaero(:,iel);

      cx = cos(xis(iel));  % xis is relative to xia.
      sx = sin(xis(iel));
      a.Ta_s(:,ic+[1:3]) = [0 0 -1;-cx sx 0;sx cx 0]; % [1 0 0;0 -cx sx;0 sx cx];

   end
end

%==========================================================
% Hub, Driveshaft, and Generator.
%==========================================================
% A direct-drive adaptation of the NREL drivetrain.  The tower-
% head mass and frequency (stiffness and inertia) and damping
% of the first drivetrain torsional mode are preserved.  Higher
% modes are not preserved. The driveshaft is made to be stiff
% in bending.

phi   = 2.5    *pi/180; % The root cone angle, for the blade coordinate
                        % system relative to the hub coordinate system.
Ld    = 5.0    /length; % Distance from rear bearing to hub center.
Lbrg  = 4.0    /length; % Bearing-to-bearing length.

Lh    = 1.5    /length; % Hub radius.
mh    = 56780  /mass;   % Hub mass.
kfact = 10;             % Stiffness factor: kfact*k(blade root) = k(hub).

Jgenr = 5.03e6 /(mass*length^2); % Additional inertia to represent the generator.

% Calibrated to L. Eliassen's Fedem model, does not perfectly
% match Jonkman.  Steel shaft of 1.5 m outer diameter and 1.429 m
% inner diameter 
% x/Ld D xis(deg) rhoA rhoJ EA EIyy EIzz GJ
drvprops = [ ...
0.0  1.5  0.0  1.282E+03  6.879E+02  3.267E+10  8.763E+09  8.763E+09  4.338E+09; ...
0.2  1.5  0.0  1.282E+03  6.879E+02  3.267E+10  8.763E+09  8.763E+09  4.338E+09; ...
0.4  1.5  0.0  1.282E+03  6.879E+02  3.267E+10  8.763E+09  8.763E+09  4.338E+09; ...
0.6  1.5  0.0  1.282E+03  6.879E+02  3.267E+10  8.763E+09  8.763E+09  4.338E+09; ...
0.8  1.5  0.0  1.282E+03  6.879E+02  3.267E+10  8.763E+09  8.763E+09  4.338E+09; ...
0.9  1.5  0.0  1.282E+03  6.879E+02  3.267E+10  8.763E+09  8.763E+09  4.338E+09; ...
1.0  1.5  0.0  1.282E+03  6.879E+02  3.267E+10  8.763E+09  8.763E+09  4.338E+09];

% Avoid numerically perfect symmetry.
drvprops(:,7) = (1 + 1e-6)*drvprops(:,7);

% Use linear interpolation.  Splines don't work well here, because
% the profile changes abruptly.
znorm      =  drvprops(:,1);
zin        = -znorm*Ld;
zn         = -[([0:Ned].')*Lbrg/Ned;Lbrg+([1:Neh].')*(Ld-Lbrg)/Neh];
ze         =  0.5*(zn(1:Ned+Neh) + zn(2:Ned+Neh+1));
Lel        =  zn(1:Ned+Neh) - zn(2:Ned+Neh+1);

sei = interp1 (zin,drvprops,ze);

Dh         = sei(:,2)/length;
xis        = sei(:,3);
rhosA      = sei(:,4)/(mass/length);
rhosJ      = sei(:,5)/(mass*length);
EA         = sei(:,6)/force;
EIyy       = sei(:,7)/(force*length^2);
EIzz       = sei(:,8)/(force*length^2);
GJ         = sei(:,9)/(force*length^2);

% Get the basic driveshaft mass (not including the hub) for later.
mdrv = sum(rhosA.*Lel);

% Additional generator mass and inertia.
igen = Ned+Neh-1;
%rhosA(igen) = rhosA(igen) + mgrot/Lel(igen);
rhosJ(igen) = rhosJ(igen) + Jgenr/Lel(igen);

% Fill in driveshaft body entries.
s.driveshaft.Ld             = Ld;
s.driveshaft.Lbrg           = Lbrg;
s.driveshaft.Lh             = Lh;
s.driveshaft.phi            = phi;
s.driveshaft.khfac          = kfact;
s.driveshaft.mh             = mh;
s.driveshaft.Jgrot          = Jgenr;

s.driveshaft.Lel(1:Ned+Neh) = Lel;
s.driveshaft.D(1:Ned+Neh)   = Dh;
s.driveshaft.xis(1:Ned+Neh) = xis;
s.driveshaft.zeta           = s.zdamp;

for iel = 1:Ned+Neh

   ic6   = 6*(iel-1);
   ic12  = 12*(iel-1);
   ic144 = 144*(iel-1);

   s.driveshaft.rhos(:,ic6+[1:6]) = ...
          diag([rhosA(iel);rhosA(iel);rhosA(iel);rhosJ(iel);0;0]);
   s.driveshaft.EEs(:,ic6+[1:6]) = ...
          diag([EA(iel);0;0;GJ(iel);EIyy(iel);EIzz(iel)]);

   s.driveshaft.ke_s(:,ic12+[1:12]) = ...
          buildkes(s.driveshaft.EEs(:,ic6+[1:6]),s.driveshaft.Lel(iel));
   s.driveshaft.me_s(:,ic12+[1:12]) = ...
          buildmes(s.driveshaft.rhos(:,ic6+[1:6]),s.driveshaft.Lel(iel));

end

for inod = 1:Ned+Neh+1
   idof = 6*(inod-1);
   s.driveshaft.Pn_B(idof+3) = zn(inod);
end

% Here are the three hub elements.
Lel     = [Lh;Lh;Lh];
xis     = [0;0;0];
rhosA   = (mh/3)./Lel;
rhosJ   = s.blade(1).rhos(4,4)*ones(3,1); % The blade root value.  Not critical.
EA      = kfact*s.blade(1).EEs(1,1)*ones(3,1);
EIyy    = kfact*s.blade(1).EEs(5,5)*ones(3,1);
EIzz    = kfact*s.blade(1).EEs(6,6)*ones(3,1);
GJ      = kfact*s.blade(1).EEs(4,4)*ones(3,1);

s.driveshaft.Lel(Ned+Neh+[1:3]) = Lel;
s.driveshaft.D(Ned+Neh+[1:3])   = a.chord(1);
s.driveshaft.xis(Ned+Neh+[1:3]) = xis;

s.driveshaft.conn = [ones(1,s.driveshaft.Nel);                              ...
             [1:s.driveshaft.Nnod-4]                                        ...
             [s.driveshaft.Nnod-3 s.driveshaft.Nnod-3 s.driveshaft.Nnod-3]; ...
             [2:s.driveshaft.Nnod-3]                                        ...
             [s.driveshaft.Nnod-2 s.driveshaft.Nnod-1 s.driveshaft.Nnod]];

for jel = 1:3

   iel   = Ned + Neh + jel;
   ic6   = 6*(iel-1);
   ic12  = 12*(iel-1);
   ic144 = 144*(iel-1);

   s.driveshaft.rhos(:,ic6+[1:6]) = ...
          diag([rhosA(jel);rhosA(jel);rhosA(jel);rhosJ(jel);0;0]);
   s.driveshaft.EEs(:,ic6+[1:6]) = ...
          diag([EA(jel);0;0;GJ(jel);EIyy(jel);EIzz(jel)]);

   s.driveshaft.ke_s(:,ic12+[1:12]) = ...
          buildkes(s.driveshaft.EEs(:,ic6+[1:6]),s.driveshaft.Lel(iel));
   s.driveshaft.me_s(:,ic12+[1:12]) = ...
          buildmes(s.driveshaft.rhos(:,ic6+[1:6]),s.driveshaft.Lel(iel));

end

idof = 6*(Ned+Neh+1);
s.driveshaft.Pn_B(idof+1)  =  Lh;
s.driveshaft.Pn_B(idof+3)  = -Ld;
s.driveshaft.Pn_B(idof+7)  =  Lh*cos(2*pi/3);
s.driveshaft.Pn_B(idof+8)  =  Lh*sin(2*pi/3);
s.driveshaft.Pn_B(idof+9)  = -Ld;
s.driveshaft.Pn_B(idof+13) =  Lh*cos(4*pi/3);
s.driveshaft.Pn_B(idof+14) =  Lh*sin(4*pi/3);
s.driveshaft.Pn_B(idof+15) = -Ld; 

%==========================================================
% Nacelle.
%==========================================================
% The masses of the nacelle elements are calibrated so that the
% total nacelle + driveshaft mass is the same as the NREL
% turbine, 240,000 kg.  This gives a "budget" from which the 
% known mass of the driveshaft is subtracted.  The remainder is
% assigned to the Nev, Ner, and Nen elements, and represents
% both the turret structure and systems.
mnac  = 240000 /mass;   % Total mass of nacelle, generator, driveshaft
                        % (not including the hub).
Lv = 1.96      /length; % Vertical distance from yaw bearing to shaft axis.
Lr = 0.2       /length; % Distance along shaft from yaw center to rear bearing.
Ln = Lbrg + Lr;         % Distance along shaft from yaw center to front bearing.

delta = 5      *pi/180; % Shaft tilt angle, rad.

% For the nacelle nose, measured from the node at the rear
% bearing.  No data available, mimic the driveshaft times a
% factor, to make it relatively stiff.  Note that the GJ
% stiffness now represents the steel cylinder, and not the
% effective stiffness of the gearbox/generator setup.
% x/Lbrg D xis(deg) rhoA rhoJ EA EIyy EIzz GJ
nacprops = [ ...
0.00  2.5  0.0  1.282E+03  6.879E+02  3.267E+10  8.763E+09  8.763E+09  1.753E+10; ...
0.25  2.5  0.0  1.282E+03  6.879E+02  3.267E+10  8.763E+09  8.763E+09  1.753E+10; ...
0.50  2.5  0.0  1.282E+03  6.879E+02  3.267E+10  8.763E+09  8.763E+09  1.753E+10; ...
0.75  2.5  0.0  1.282E+03  6.879E+02  3.267E+10  8.763E+09  8.763E+09  1.753E+10; ...
1.00  2.5  0.0  1.282E+03  6.879E+02  3.267E+10  8.763E+09  8.763E+09  1.753E+10];

NevNen     =  Nev + Ner + Nen;
znorm      =  nacprops(:,1);
zin        =  znorm*(Lv + Ln);
zn         =  [([0:Nev].')*Lv/Nev;      ...
               Lv + ([1:Ner].')*Lr/Ner; ...
               Lv + Lr + ([1:Nen].')*(Ln - Lr)/Nen];
ze         =  0.5*(zn(1:NevNen) + zn(2:NevNen+1));
Lel        =  -zn(1:NevNen) + zn(2:NevNen+1);

sei = interp1 (zin,nacprops,ze);

Dh         = sei(:,2)/length;
xis        = sei(:,3);
rhosA      = sei(:,4)/(mass/length);
rhosJ      = sei(:,5)/(mass*length);
EA         = sei(:,6)/force;
EIyy       = sei(:,7)/(force*length^2);
EIzz       = sei(:,8)/(force*length^2);
GJ         = sei(:,9)/(force*length^2);

% Assign remaining mass to the nacelle structure.
mrem  = mnac - mdrv ...
      - sum(rhosA(Nev+Ner+1:NevNen).*Lel(Nev+Ner+1:NevNen));

rhosA(1:Nev+Ner) = mrem/sum(Lel(1:Nev+Ner)); 

s.nacelle.Lv     = Lv;
s.nacelle.Lr     = Lr;
s.nacelle.Ln     = Ln;
s.nacelle.delta  = delta;
s.nacelle.D(:)   = Dh;

s.nacelle.Lel(:) = Lel;
s.nacelle.xis(:) = xis;
s.nacelle.zeta   = s.zdamp;

s.nacelle.conn = [ones(1,s.nacelle.Nel); ...
                  [1:s.nacelle.Nel];     ...
                  [2:s.nacelle.Nel+1]];

for iel = 1:Nev+Ner+Nen

   ic6  = 6*(iel-1);
   ic12 = 12*(iel-1);
   ic144 = 144*(iel-1);

   s.nacelle.rhos(:,ic6+[1:6]) = ...
          diag([rhosA(iel);rhosA(iel);rhosA(iel);rhosJ(iel);0;0]);
   s.nacelle.EEs(:,ic6+[1:6]) = ...
          diag([EA(iel);0;0;GJ(iel);EIyy(iel);EIzz(iel)]);

   s.nacelle.ke_s(:,ic12+[1:12]) = ...
          buildkes(s.nacelle.EEs(:,ic6+[1:6]),s.nacelle.Lel(iel));
   s.nacelle.me_s(:,ic12+[1:12]) = ...
          buildmes(s.nacelle.rhos(:,ic6+[1:6]),s.nacelle.Lel(iel));

end

for inod = 1:Nev+1
   idof = 6*(inod-1);
   s.nacelle.Pn_B(idof+3) = zn(inod);
end
for inod = Nev+1:Nev+Ner+Nen+1
   idof = 6*(inod-1);
   s.nacelle.Pn_B(idof+1) = -(zn(inod) - zn(Nev+1))*cos(delta);
   s.nacelle.Pn_B(idof+3) =  zn(Nev+1) + (zn(inod) - zn(Nev+1))*sin(delta);
end

%==========================================================
% Tower.
%==========================================================

% Split the tower between the "tower" and "foundation" in STAS.
Lt = 43.8      /length; % 38.8/length; % 
Cdt0 = 0.80;            % Over the nominal 0.65, margin to account for icing
                        % and weathering of the surface.  Ref DNV-OS-J101.

% Tower.
% x/Lt Dt xis(deg) rhoA rhoJ EA EIyy EIzz GJ
towprops = [ ...
0.0 4.94 0 3916 2.20e4 9.68e10 2.91e11 2.91e11 2.24e11; ...
0.2 4.72 0 3616 1.85e4 8.94e10 2.46e11 2.46e11 1.89e11; ...
0.4 4.51 0 3329 1.56e4 8.23e10 2.06e11 2.06e11 1.59e11; ...
0.6 4.30 0 3053 1.30e4 7.54e10 1.72e11 1.72e11 1.32e11; ...
0.8 4.08 0 2789 1.07e4 6.89e10 1.42e11 1.42e11 1.09e11; ...
1.0 3.87 0 2536 8.74e3 6.27e10 1.16e11 1.16e11 8.91e10];

% Avoid numerically perfect symmetry.
towprops(:,7) = (1 + 1e-6)*towprops(:,7);

% Use linear interpolation.  Splines don't work well here, because
% the profile changes abruptly.
znorm      =  towprops(:,1);
zin        =  znorm*Lt;
zn         =  ([0:Net].')*Lt/Net;
ze         =  0.5*(zn(1:Net) + zn(2:Net+1));
Lel        = -zn(1:Net) + zn(2:Net+1);

sei = interp1 (zin,towprops,ze);

Dh         = sei(:,2)/length;
xis        = sei(:,3);
rhosA      = sei(:,4)/(mass/length);
rhosJ      = sei(:,5)/(mass*length);
EA         = sei(:,6)/force;
EIyy       = sei(:,7)/(force*length^2);
EIzz       = sei(:,8)/(force*length^2);
GJ         = sei(:,9)/(force*length^2);

% Fill in tower body entries.
s.tower.Lt             = Lt;
s.tower.Cd(:)          = Cdt0;

s.tower.Lel(:)         = Lel;
s.tower.D(:)           = Dh;
s.tower.xis(:)         = xis;
s.tower.zeta           = s.zdamp;

s.tower.conn = [ones(1,s.tower.Nel); ...
                [1:s.tower.Nel];     ...
                [2:s.tower.Nel+1]];

for iel = 1:Net

   ic6  = 6*(iel-1);
   ic12 = 12*(iel-1);
   ic144 = 144*(iel-1);

   s.tower.rhos(:,ic6+[1:6]) = ...
          diag([rhosA(iel);rhosA(iel);rhosA(iel);rhosJ(iel);0;0]);
   s.tower.EEs(:,ic6+[1:6]) = ...
          diag([EA(iel);0;0;GJ(iel);EIyy(iel);EIzz(iel)]);

   s.tower.ke_s(:,ic12+[1:12]) = ...
          buildkes(s.tower.EEs(:,ic6+[1:6]),s.tower.Lel(iel));
   s.tower.me_s(:,ic12+[1:12]) = ...
          buildmes(s.tower.rhos(:,ic6+[1:6]),s.tower.Lel(iel));

end

for inod = 1:Net+1
   idof = 6*(inod-1);
   s.tower.Pn_B(idof+3) = zn(inod);
end

%==========================================================
% Foundation.
%==========================================================
% The foundation consists arbitrarily of the lower half of
% the tower.
Lf = 43.8      /length; % 38.8/length; % 

rhow = 1024  /(mass/length^3);
Cdf0 = 1.0;            % To account for marine growth.
wel = [6:20]';         % Must match length of Nwel at top of input file.

% x/Lt Dt xis(deg) rhoA rhoJ EA EIyy EIzz GJ
fnprops = [ ...
0.0 6.00 0 5591 4.64e4 1.38e11 6.14e11 6.14e11 4.73e11; ...
0.2 5.79 0 5232 4.03e4 1.29e11 5.35e11 5.35e11 4.11e11; ...
0.4 5.57 0 4886 3.50e4 1.21e11 4.63e11 4.63e11 3.57e11; ...
0.6 5.36 0 4551 3.01e4 1.12e11 3.99e11 3.99e11 3.07e11; ...
0.8 5.15 0 4228 2.58e4 1.04e11 3.42e11 3.42e11 2.63e11; ...
1.0 4.94 0 3916 2.20e4 9.68e10 2.91e11 2.91e11 2.24e11];

% Avoid numerically perfect symmetry.
fnprops(:,7) = (1 + 1e-6)*fnprops(:,7);

% Use linear interpolation.
znorm      =  fnprops(:,1);
zin        =  znorm*Lf;
zn         =  ([0:Nef].')*Lf/Nef;
ze         =  0.5*(zn(1:Nef) + zn(2:Nef+1));
Lel        =  -zn(1:Nef) + zn(2:Nef+1);

sei = interp1 (zin,fnprops,ze);

Dh         = sei(:,2)/length;
xis        = sei(:,3);
rhosA      = sei(:,4)/(mass/length);
rhosJ      = sei(:,5)/(mass*length);
EA         = sei(:,6)/force;
EIyy       = sei(:,7)/(force*length^2);
EIzz       = sei(:,8)/(force*length^2);
GJ         = sei(:,9)/(force*length^2);

% Soil spring properties.



%                                  --------------------
scale = 0.0;  % <================= | No soil springs. |
%                                  --------------------



dkxy = scale*[1.15e8; 9.04e7; 6.94e7; 5.12e7; 3.25e7; 1.89e7]/(force/length);

dkz = dkxy*(1.62/2.06);

kxg = zeros(Nef+1,1);
kyg = zeros(Nef+1,1);
kzg = zeros(Nef+1,1);
cxg = zeros(Nef+1,1);
cyg = zeros(Nef+1,1);
czg = zeros(Nef+1,1);
kthzg = zeros(Nef+1,1);
cthzg = zeros(Nef+1,1);

% Fill in foundation body entries.
s.foundation.Lf             = Lf;
s.foundation.wel            = wel;
s.foundation.Nwater         = Nwater;
s.foundation.rhow           = rhow;
s.foundation.Cd(:)          = Cdf0;
s.foundation.k(:,:)         = [kxg kyg kzg kthzg cxg cyg czg cthzg].';

s.foundation.Lel(:)         = Lel;
s.foundation.D(:)           = Dh;
s.foundation.xis(:)         = xis;
s.foundation.zeta           = s.zdamp;

s.foundation.conn = [ones(1,s.foundation.Nel); ...
                     [1:s.foundation.Nel];     ...
                     [2:s.foundation.Nel+1]];

for iel = 1:Nef

   ic6   = 6*(iel-1);
   ic12  = 12*(iel-1);
   ic144 = 144*(iel-1);

   s.foundation.rhos(:,ic6+[1:6]) = ...
          diag([rhosA(iel);rhosA(iel);rhosA(iel);rhosJ(iel);0;0]);
   s.foundation.EEs(:,ic6+[1:6]) = ...
          diag([EA(iel);0;0;GJ(iel);EIyy(iel);EIzz(iel)]);

   s.foundation.ke_s(:,ic12+[1:12]) = ...
          buildkes(s.foundation.EEs(:,ic6+[1:6]),s.foundation.Lel(iel));
   s.foundation.me_s(:,ic12+[1:12]) = ...
          buildmes(s.foundation.rhos(:,ic6+[1:6]),s.foundation.Lel(iel));

end

for inod = 1:Nef+1
   idof = 6*(inod-1);
   s.foundation.Pn_B(idof+3) = zn(inod) - zn(1);
end

