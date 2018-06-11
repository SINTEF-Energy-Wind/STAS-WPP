function [s,a] = STASTurbine_DTU10MW ()

% -------------------------------------
% 10 MW Offshore Reference Wind Turbine
% -------------------------------------
%
% This function defines input properties for the DTU 10 MW Reference
% Wind Turbine.  The nacelle and drivetrain properties are adapted
% from the NOWITECH 10 MW reference wind turbine.
%
% Version:        Changes:
% --------        -------------
% 31.08.2017      Original code.  A complete rewrite of the original
%                 STASTurbine.m in order to include data structures.
%
% Version:        Verification:
% --------        -------------
% 31.08.2017      
%
% Properties:
% -----------
%

%==========================================================
% Normalization.
%==========================================================
length = 1
time   = 1
power  = 1e6;
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
Nef    = 20;
Nmud   = 5;  % (Nmud is part of Nef, specifies number below seabed.)
Nwater = 10; % (Nwater  "      "   , number between seabed and waterline.)
Net    = 10;
Nev    = 2;
Ner    = 2;
Nen    = 4;
Ned    = 4;
Neh    = 2;
Neb    = 16;

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

s.zdamp = 0; % 0.008; % Modal structural damping ratio.
                      % or ...
s.adamp = 0;          % Mass-proportional damping factor.
s.bdamp = 0.0002;     % Stiffness-proportional damping factor.

% Correct zero-lift aoa values for inboard DTU 10 MW profiles.
load('-binary','aoazs_DTU10MW.bin');
Nfoils = size(aoazs,1);
aoazs(1) = 0;
aoazs(2) = 4.73*pi/180;
aoazs(3) = -2.50*pi/180;
save('-binary','aoazs_DTU10MW.bin','aoazs');
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
load('DTU10MWAeroProfile.txt');
load('DTU10MWBladeStructure.txt');
Nain = size(DTU10MWAeroProfile,1);      % Number of rows in the input files.
Nsin = size(DTU10MWBladeStructure,1);

Lb   = 86.466   /length; % Blade length from pitch bearing flange to tip.

% Get the radial coordinates at which the aero and structural data
% is given.
rbain = Lb*DTU10MWAeroProfile(:,1)/length;
rbsin = Lb*DTU10MWBladeStructure(:,1)/length;

% Distribute nodes along the blades.  Refined near the root and tip,
% coarser in the middle.
rbnod   = distributer (Neb+1,rbain(1),rbain(Nain),-0.35*pi,0.45*pi);
rbel    = 0.5*(rbnod(2:Neb+1) + rbnod(1:Neb));

% Compute spline matrices for later interpolation.  If we have some
% quantity gi at input table points, then spline coefficients are
% computed by kgi = AB*gi, and interpolated values by g = S*kgi.
[Aaero,Baero]     = splineCoefficients (rbain);
ABaero            = Aaero\Baero;
[Saero,dSaero]    = splineSMatrix (rbain,rbel);
[Snaero,dSnaero]  = splineSMatrix (rbain,rbnod);
[Astruct,Bstruct] = splineCoefficients (rbsin);
ABstruct          = Astruct\Bstruct;
[Sstruct,dSstruct]= splineSMatrix (rbsin,rbel);

% Structural properties.  The input tables from DTU are given in the
% hub coordinate system, that is, zero root cone.
%
% In the structural calculations the structural twist defines the
% section coordinate system, while in the aerodynamic calculations
% the aerodynamic twist defines the airfoil coordinate system.
%
% This input is simplified, for equivalent beam properties.  STAS
% accepts the full 6-by-6 inertia and stiffness section property
% matrices as input, via rhos and EEs, if these properties are
% available.  At present I've ignored the elastic center and shear
% center, but have included the center-of-mass offset in the chordwise
% direction.
%
% Spline coefficients of structural properties, and interpolated 
% values at the desired element locations.
ksei    = ABstruct*DTU10MWBladeStructure;
sei     = Sstruct*ksei;

EA      =  sei(:,4)/force;
EIyy    =  sei(:,5)/(force*length^2);
EIzz    =  sei(:,6)/(force*length^2);
GJ      =  sei(:,7)/(force*length^2);
egy     = -sei(:,8)/length;       % Centroid in structural section CS.
egz     =  sei(:,9)/length;       % My y^s=-x_in, z^s=y_in.
xis     = -sei(:,10)*pi/180;      % My structural twist uses the opposite CS.
rhosA   =  (sei(:,13) + sei(:,19))/(mass/length); % Mass per unit length.
rhosJ   =  sei(:,14)/(mass*length); % Torsional inertia per unit length.
cgy     = -sei(:,17)/length;      % Center of mass in structural section CS.
cgz     =  sei(:,18)/length;

% Aerodynamic properties.
kai     = ABaero*DTU10MWAeroProfile;
ai      = Saero*kai;             % Values at elements.
ani     = Snaero*kai;            % Values at nodes.

chord   = ai(:,2)/(1000*length);
xia     = ai(:,3)*pi/180;        % Aero twist, rad.  Positive is LE into wind.
xpc     = ai(:,4);               % Origin of section CS aft of the LE,
                                 % fraction of chord.
tca     = ai(:,5)/100;           % t/c.

% Element coordinates, pitch CS.
xbaero  = [rbel ai(:,7)/(1000*length) ai(:,6)/(1000*length)].';

acent   = ai(:,8)/100;           % Aero center aft of LE, fraction of chord.

% Nodal coordinates, pitch CS.
xbnaero = [rbnod ani(:,7)/(1000*length) ani(:,6)/(1000*length)].'; 
Lel     = norm(xbnaero(:,2:Neb+1) - xbnaero(:,1:Neb),2,'cols').';

% Compute foil weights from t/c.
load('-ascii','tcfoils_DTU10MW.txt');
foilwt = getFoilWeights (tcfoils_DTU10MW,tca);

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

load('-binary','aoas_DTU10MW.bin');
load('-binary','kfoils_DTU10MW.bin');
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
% The drivetrain is based on the NOWITECH Reference Turbine,
% whose drivetrain definition is more complete than that of
% the DTU turbine.
%
% The hub mass is that of the DTU Reference Turbine.

phi  = 2.5     *pi/180; % The root cone angle, for the blade coordinate
                        % system relative to the hub coordinate system.
Ld    = 8.0    /length; % Distance from rear bearing to hub center,
                        % = 5.4 m driveshaft to flange + about 2.6 m
                        % to hub center.
Lbrg  = 3.714  /length; % Bearing-to-bearing length.

Lh    = 2.8    /length; % Hub radius.
mh    = 105520 /mass;   % Hub mass.
kfact = 5;              % Stiffness factor: kfact*k(blade root) = k(hub).

rrot  = 5.970  /length; % Effective r of generator rotor for Jgen calculation.
mgen  = 240000 /mass;   % Memo AN 15.12.35, SINTEF Energi, 2015.
mgrot = 0.5*mgen;       % Generator rotor mass.
Jgenr = 0.625*mgrot*(rrot^2); % Rough approximation: Half in housing/  


% The driveshaft.
% The final element is located inside the hub, and its stiffness
% and mass properties should reflect this.  Mass: hub mass is
% accounted for by the three radial elements, so let the last
% driveshaft element be light.  Stiffness: should be stiff.
% x/Ld D xis(deg) rhoA rhoJ EA EIyy EIzz GJ
drvprops = [ ...
0.000000E+00  1.700000E+00  0.000000E+00  5.952729E+03  3.506529E+03  1.425301E+11  4.197958E+10  4.197958E+10  3.308978E+10; ...
5.177557E-02  1.700000E+00  0.000000E+00  5.952729E+03  3.506529E+03  1.425301E+11  4.197958E+10  4.197958E+10  3.308978E+10; ...
2.447354E-01  1.700000E+00  0.000000E+00  5.952729E+03  3.506529E+03  1.425301E+11  4.197958E+10  4.197958E+10  3.308978E+10; ...
2.857820E-01  1.760000E+00  0.000000E+00  6.653673E+03  4.160209E+03  1.593133E+11  4.980532E+10  4.980532E+10  3.925831E+10; ...
3.268285E-01  1.820000E+00  0.000000E+00  7.384730E+03  4.892845E+03  1.768175E+11  5.857631E+10  5.857631E+10  4.617192E+10; ...
3.678750E-01  1.880000E+00  0.000000E+00  8.145898E+03  5.710275E+03  1.950426E+11  6.836244E+10  6.836244E+10  5.388569E+10; ...
4.089215E-01  1.940000E+00  0.000000E+00  8.937179E+03  6.618540E+03  2.139888E+11  7.923604E+10  7.923604E+10  6.245664E+10; ...
4.499680E-01  2.000000E+00  0.000000E+00  9.758572E+03  7.623885E+03  2.336560E+11  9.127186E+10  9.127186E+10  7.194370E+10; ...
5.302202E-01  2.000000E+00  0.000000E+00  9.758572E+03  7.623885E+03  2.336560E+11  9.127186E+10  9.127186E+10  7.194370E+10; ...
5.562038E-01  2.000000E+00  0.000000E+00  9.758572E+03  7.623885E+03  2.336560E+11  9.127186E+10  9.127186E+10  7.194370E+10; ...
5.749964E-01  2.000000E+00  0.000000E+00  1.427540E+04  9.707270E+03  3.418053E+11  1.162138E+11  1.162138E+11  9.160382E+10; ...
6.174716E-01  2.000000E+00  0.000000E+00  1.427540E+04  9.707270E+03  3.418053E+11  1.162138E+11  1.162138E+11  9.160382E+10; ...
6.395241E-01  2.000000E+00  0.000000E+00  1.427540E+03  9.707270E+02  3.418053E+12  1.162138E+12  1.162138E+12  9.160382E+11; ...
6.750000E-01  2.000000E+00  0.000000E+00  1.427540E+03  9.707270E+02  3.418053E+12  1.162138E+12  1.162138E+12  9.160382E+11; ...
1.000000E+00  2.000000E+00  0.000000E+00  1.427540E+03  9.707270E+02  3.418053E+12  1.162138E+12  1.162138E+12  9.160382E+11];

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

% Additional generator mass and inertia.  The generator attaches
% to the driveshaft at the node outboard of the front bearing,
% just behind the hub.  However, for studies, let this be selected
% by genFlag, from STASElectric.m.
igen = Ned+Neh-1;
rhosA(igen) = rhosA(igen) + mgrot/Lel(igen);
rhosJ(igen) = rhosJ(igen) + Jgenr/Lel(igen);

% Fill in driveshaft body entries.
s.driveshaft.Ld             = Ld;
s.driveshaft.Lbrg           = Lbrg;
s.driveshaft.Lh             = Lh;
s.driveshaft.phi            = phi;
s.driveshaft.khfac          = kfact;
s.driveshaft.mh             = mh;
s.driveshaft.rrot           = rrot;
s.driveshaft.mgrot          = mgrot;
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
% The masses of the nacelle elements are calibrated so that
% the total nacelle + driveshaft mass is the same as the
% DTU turbine, 446,026 kg.  This gives a "budget" from which
% the known masses of the driveshaft and nacelle nose are
% subtracted.  The remainder is assigned to the Nev and Ner
% elements, and represents both the turret structure and
% systems.
mnac  = 446026 /mass;   % Total mass of nacelle, generator, driveshaft
                        % (not including the hub).
Lv    = 2.75   /length; % Vertical distance from yaw bearing to shaft axis.
Lr    = 2      /length; % Distance along shaft from yaw center to rear bearing.
Ln    = Lbrg + Lr;      % Distance along shaft from yaw center to front bearing.

delta = 5      *pi/180; % Shaft tilt angle, rad.
rsta  = 5.897  /length; % Effective r of generator stator for Jgen calculation.
mgsta = 0.5*mgen;       % Generator stator mass.
Jgens = 0.625*mgsta*(rsta^2); % support at r/2, half in outer part at r.

% The nacelle turret and nose are parameterized beginning at the 
% yaw bearing, with elements arranged vertically up to the 
% intersection with the driveshaft axis, then horizontally (with
% tilt delta) to the nose.  Nev elements are vertical, Ner to
% the rear bearing (nose flange), and Nen to the front bearing.
% x/Lbrg D xis(deg) rhoA rhoJ EA EIyy EIzz GJ
nacprops = [ ...
0.0000E+00  5.00  0.0000E+00  0.0000E+00  4.0000E+04  4.0000E+11  4.0000E+11  4.0000E+11  8.0000E+11; ...
5.0000E-01  5.00  0.0000E+00  0.0000E+00  4.0000E+04  4.0000E+11  4.0000E+11  4.0000E+11  8.0000E+11; ...
5.6120E-01  3.13  0.0000E+00  1.3375E+04  2.8748E+04  3.2024E+11  3.4417E+11  3.4417E+11  6.8833E+11; ...
6.7090E-01  3.13  0.0000E+00  1.3375E+04  2.8748E+04  3.2024E+11  3.4417E+11  3.4417E+11  6.8833E+11; ...
7.8060E-01  3.13  0.0000E+00  6.7585E+03  1.5529E+04  1.6182E+11  1.8591E+11  1.8591E+11  3.7183E+11; ...
8.9030E-01  3.13  0.0000E+00  1.9078E+04  3.8568E+04  4.5680E+11  4.6173E+11  4.6173E+11  9.2345E+11; ...
1.0000E+00  3.13  0.0000E+00  1.0286E+04  2.2820E+04  2.4627E+11  2.7320E+11  2.7320E+11  5.4640E+11];

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

% Add the generator stator mass to the outermost nose element.
igen = Nev + Ner + Nen;
rhosA(igen) = rhosA(igen) + mgsta/Lel(igen);
rhosJ(igen) = rhosJ(igen) + Jgens/Lel(igen);

% Assign remaining mass to the nacelle structure.
mrem  = mnac - mgrot - mdrv ...
      - sum(rhosA(Nev+Ner+1:NevNen).*Lel(Nev+Ner+1:NevNen));
rhosA(1:Nev+Ner) = mrem/sum(Lel(1:Nev+Ner)); 

s.nacelle.Lv     = Lv;
s.nacelle.Lr     = Lr;
s.nacelle.Ln     = Ln;
s.nacelle.delta  = delta;
s.nacelle.D(:)   = Dh;
s.nacelle.rstat  = rsta;
s.nacelle.mgstat = mgsta;
s.nacelle.Jgstat = Jgens;

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
% The tower, starting at 10 m elevation.
Lt   = 105.63  /length; % Originally 115.63.  Compressed length by 10 m
                        % to accommodate the above-waterline portion of
                        % the piled foundation.
Cdt0 = 0.80;            % Over the nominal 0.65, margin to account for icing
                        % and weathering of the surface.  Ref DNV-OS-J101.

% x/Lt Dt xis(deg) rhoA rhoJ EA EIyy EIzz GJ
towprops = [ ...
0.0000E+00  8.3000E+00  0.0000E+00  1.4207E+04  2.4060E+05  3.8007E+11  3.2182E+12  3.2182E+12  2.4826E+12; ...
9.9455E-02  8.0215E+00  0.0000E+00  1.4207E+04  2.4060E+05  3.8007E+11  3.2182E+12  3.2182E+12  2.4826E+12; ...
1.9891E-01  7.7430E+00  0.0000E+00  1.2791E+04  2.0236E+05  3.4217E+11  2.7067E+12  2.7067E+12  2.0880E+12; ...
2.9837E-01  7.4646E+00  0.0000E+00  1.1439E+04  1.6867E+05  3.0602E+11  2.2561E+12  2.2561E+12  1.7404E+12; ...
3.9782E-01  7.1861E+00  0.0000E+00  1.0153E+04  1.3916E+05  2.7160E+11  1.8614E+12  1.8614E+12  1.4360E+12; ...
4.9728E-01  6.9076E+00  0.0000E+00  8.9315E+03  1.1349E+05  2.3893E+11  1.5181E+12  1.5181E+12  1.1711E+12; ...
5.9673E-01  6.6292E+00  0.0000E+00  7.7754E+03  9.1321E+04  2.0800E+11  1.2215E+12  1.2215E+12  9.4229E+11; ...
6.9619E-01  6.3507E+00  0.0000E+00  6.6844E+03  7.2329E+04  1.7882E+11  9.6746E+11  9.6746E+11  7.4633E+11; ...
7.9564E-01  6.0722E+00  0.0000E+00  5.6586E+03  5.6214E+04  1.5138E+11  7.5190E+11  7.5190E+11  5.8004E+11; ...
8.9510E-01  5.7937E+00  0.0000E+00  4.6980E+03  4.2684E+04  1.2568E+11  5.7093E+11  5.7093E+11  4.4043E+11; ...
1.0000E+00  5.5000E+00  0.0000E+00  3.8025E+03  3.1465E+04  1.0172E+11  4.2087E+11  4.2087E+11  3.2467E+11];

% Avoid numerically perfect symmetry.
towprops(:,7) = (1 + 1e-6)*towprops(:,7);

% Use linear interpolation.  Splines don't work well here, because
% the profile changes abruptly.
znorm      =  towprops(:,1);
zin        =  znorm*Lt;
zn         =  ([0:Net].')*Lt/Net;
ze         =  0.5*(zn(1:Net) + zn(2:Net+1));
Lel        = -zn(1:Net) + zn(2:Net+1);

sei = interp1 (zn,towprops,ze);

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
% The foundation consists of a monopile driven 42 m into the
% seabed, and extending to an elevation of 20 m above the
% seabed (and 10 m below sealevel); and a 20 m high transition
% piece, from -10 m to +10 m relative to sealevel.
Lf = 82      /length;
Lt = 105.63  /length;  % Originally 115.63.  Compressed length by 10 m
                       % to accommodate the above-waterline portion of
                       % the piled foundation.

rhow = 1024  /(mass/length^3);
Cdf0 = 1.0;            % To account for marine growth.
wel = [6:20]';         % Must match length of Nwel at top of input file.

% Based on a design obtained from L. Eliassen, using soil P-Y
% relationships representative of the Dogger Bank.  Do = 9 m.
% Includes refinement near the waterline, for wave loads.  The
% change in properties corresponds to the change from the
% monopile to transition piece.
% x/Lf Df xis(deg) rhoA rhoJ EA EIyy EIzz GJ
fnprops = [ ...
0.0000E+00  9.0000E+00  0.0000E+00  2.1949E+04  4.3469E+05  5.8716E+11  5.8144E+12  5.8144E+12  4.4854E+12; ...
1.0244E-01  9.0000E+00  0.0000E+00  2.1949E+04  4.3469E+05  5.8716E+11  5.8144E+12  5.8144E+12  4.4854E+12; ...
2.0488E-01  9.0000E+00  0.0000E+00  2.1949E+04  4.3469E+05  5.8716E+11  5.8144E+12  5.8144E+12  4.4854E+12; ...
3.0732E-01  9.0000E+00  0.0000E+00  2.1949E+04  4.3469E+05  5.8716E+11  5.8144E+12  5.8144E+12  4.4854E+12; ...
4.0976E-01  9.0000E+00  0.0000E+00  2.1949E+04  4.3469E+05  5.8716E+11  5.8144E+12  5.8144E+12  4.4854E+12; ...
5.1220E-01  9.0000E+00  0.0000E+00  2.1949E+04  4.3469E+05  5.8716E+11  5.8144E+12  5.8144E+12  4.4854E+12; ...
5.6098E-01  9.0000E+00  0.0000E+00  2.1949E+04  4.3469E+05  5.8716E+11  5.8144E+12  5.8144E+12  4.4854E+12; ...
6.0976E-01  9.0000E+00  0.0000E+00  2.1949E+04  4.3469E+05  5.8716E+11  5.8144E+12  5.8144E+12  4.4854E+12; ...
6.5854E-01  9.0000E+00  0.0000E+00  2.1949E+04  4.3469E+05  5.8716E+11  5.8144E+12  5.8144E+12  4.4854E+12; ...
7.0732E-01  9.0000E+00  0.0000E+00  2.1949E+04  4.3469E+05  5.8716E+11  5.8144E+12  5.8144E+12  4.4854E+12; ...
7.5610E-01  9.0000E+00  0.0000E+00  3.4881E+04  6.8168E+05  9.3313E+11  9.1180E+12  9.1180E+12  7.0339E+12; ...
7.8049E-01  9.0000E+00  0.0000E+00  3.4743E+04  6.7361E+05  9.2943E+11  9.0101E+12  9.0101E+12  6.9507E+12; ...
8.0488E-01  9.0000E+00  0.0000E+00  3.4467E+04  6.5768E+05  9.2204E+11  8.7970E+12  8.7970E+12  6.7862E+12; ...
8.2927E-01  9.0000E+00  0.0000E+00  3.4191E+04  6.4200E+05  9.1466E+11  8.5872E+12  8.5872E+12  6.6244E+12; ...
8.5366E-01  9.0000E+00  0.0000E+00  3.3914E+04  6.2657E+05  9.0727E+11  8.3808E+12  8.3808E+12  6.4652E+12; ...
8.7805E-01  9.0000E+00  0.0000E+00  3.3638E+04  6.1139E+05  8.9988E+11  8.1778E+12  8.1778E+12  6.3086E+12; ...
9.0244E-01  9.0000E+00  0.0000E+00  3.3362E+04  5.9645E+05  8.9249E+11  7.9780E+12  7.9780E+12  6.1545E+12; ...
9.2683E-01  9.0000E+00  0.0000E+00  3.3086E+04  5.8176E+05  8.8510E+11  7.7816E+12  7.7816E+12  6.0029E+12; ...
9.5122E-01  9.0000E+00  0.0000E+00  3.2810E+04  5.6732E+05  8.7771E+11  7.5883E+12  7.5883E+12  5.8539E+12; ...
9.7561E-01  9.0000E+00  0.0000E+00  3.2533E+04  5.5311E+05  8.7032E+11  7.3983E+12  7.3983E+12  5.7073E+12; ...
1.0000E+00  9.0000E+00  0.0000E+00  3.2257E+04  5.3915E+05  8.6293E+11  7.2116E+12  7.2116E+12  5.5632E+12];

% Avoid numerically perfect symmetry.
fnprops(:,7) = (1 + 1e-6)*fnprops(:,7);

% Use an exact number of elements here, in order to get the correct
% resolution near the waterline etc.
znorm      =  fnprops(:,1);
zn        =  -Lf + znorm*Lf;
ze         =  0.5*(zn(1:Nef) + zn(2:Nef+1));
Lel        =  -zn(1:Nef) + zn(2:Nef+1);

sei = interp1 (zn,fnprops,ze);

Dh         = sei(:,2)/length;
xis        = sei(:,3);
rhosA      = sei(:,4)/(mass/length);
rhosJ      = sei(:,5)/(mass*length);
EA         = sei(:,6)/force;
EIyy       = sei(:,7)/(force*length^2);
EIzz       = sei(:,8)/(force*length^2);
GJ         = sei(:,9)/(force*length^2);

CmA = [0.50*pi*(Dh(Nmud+[1:Nwater]).^2); ... % Ca = 1, + flooded tower below waterline.
       0.25*pi*(Dh(Nmud+Nwater+[1:(Nwel-Nwater)]).^2)]; 

%CmA = 0.25*pi*(Dh(Nmud+1:Nmud+Nwel).^2); % Ca = 1.

% Increase the effective rho*A density to account for soil
% within the pile.  Density of 2000 kg/m^3.
rhosA(1:Nmud) = rhosA(1:Nmud) + 0.25*pi*(Dh(1:Nmud).^2) ...
              * 2e3/(mass/length^3);

% Soil spring properties.
% dkxy is based on values obtained from L. Eliassen,
% representative of the Dogger Bank.  dkz is estimated
% based on previous calculations (E. Van Buren's code)
% on dense sand and clay.  It is scaled from the given
% values of dkxy.
% depth = -42  -33.6  -25.2   -16.8  -8.4   0
scale = 1.0;
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

% Apply stiffness to the bottom nodes, corresponding to
% the bottom Nmud (submerged in the seabed) elements.
kxg(1) = 0.5*dkxy(1)*Lel(1);
kxg(2:Nmud) = 0.5*dkxy(2:Nmud) ...
           .* (Lel(1:Nmud-1) + Lel(2:Nmud));
kxg(Nmud+1) = 0.5*dkxy(Nmud+1)*Lel(Nmud);
kyg = kxg;
kzg(1:Nmud+1) = kxg(1:Nmud+1).*dkz./dkxy;

% Damgaard 2013, should be a bit under 0.01 damping ratio for
% the first mode.  If f = 0.25 Hz, omega = 1.57 rad/s, then 
% c = 2 zeta k/omega.  Say zeta = 0.008, then c = 0.01 k.
% This isn't rigorously correct, and it should be verified
% what the actual modal damping is.
dfac = 0.01;
cxg = dfac*kxg;
cyg = dfac*kyg;
czg = dfac*kzg;

% A torsional stiffness and damping are needed as well.
% These can be derived as R^2 times the tension stiffness.
kthzg = [(Dh.^2).*kzg(1:Nef)/4; (Dh(Nef)^2)*kzg(Nef+1)/4];
cthzg = [(Dh.^2).*czg(1:Nef)/4; (Dh(Nef)^2)*czg(Nef+1)/4];

% Fill in foundation body entries.
s.foundation.Lf             = Lf;
s.foundation.wel            = wel;
s.foundation.Nwater         = Nwater;
s.foundation.rhow           = rhow;
s.foundation.CmA(:)         = CmA;
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



