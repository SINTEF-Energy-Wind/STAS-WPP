function [s,a] = STASTurbine_Tjaereborg ()

% -------------------------------------
% Tjæreborg 2 MW prototype wind turbine
% -------------------------------------
%
% 
%
% Version:        Changes:
% --------        -------------
% 01.03.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 01.03.2018      
%
% Properties:
% -----------
%

%==========================================================
% Normalization.
%==========================================================
length = 1;
time   = 1;
power  = 1e6; % 1; % 
voltage = sqrt(power);
%---------------------------------------
velocity = length/time;
mass   = power*(time^3)/(length^2);
force  = mass*length/(time^2);
stress = force/(length^2);
ndens  = mass/(length^3);
nvisc  = mass/(length*time);
stiffness = force/length;
damping = force*time/length;
current = power/voltage;
resistance = voltage/current;
inductance = voltage*time/current;
capacitance = current*time/voltage;
flux   = voltage*time;

LTMnorms = [length;time;mass;current];
save('-ascii','LTMnorms.txt','LTMnorms');

%==========================================================
% General.
%==========================================================
Nef    = 4;
Net    = 4;
Nev    = 2;
Ner    = 2;
Nen    = 4;
Ned    = 4;
Neh    = 2;
Neb    = 12;

% Unused here...
Nmud   = Nef/2;  % (Nmud is part of Nef, specifies number below seabed.)
Nwater = Nef/2;  % (Nwater  "      "   , number between seabed and waterline.)

% Number of modes to retain.  Retain as few as possible in order
% that the numerical conditioning of the structural matrices
% will be as good as possible.  (This can also be accomplished
% by later partitioning of the system matrices, if, for instance,
% a study is focused on a group of high-frequency modes.)
%              |  Caution, if you use more than this number, the
%              V  code has not been validated.
Nfnd = 4;  %   20;  % Foundation.
Ntow = 4;  %   20;  % Tower.
Nnac = 2;  %   8;   % Nacelle.
Ndrv = 2;  %   18;  % Driveshaft.
Nbld = 16; %   20;  % Blades.

Nwel = 0;

% Set to positive integers if icp contains the control points
% for a spline representation of the aerodynamic states. Set
% to negative integers if icp contains the blade body modes
% for a modal-slave representation of the aerodynamic states.
% In the latter case, the selection of modes has to be verified
% as the turbine design is altered.
%icp = round(1 + (Neb-1)*[0.00 0.25 0.50 0.70 0.85 1.00]'); % Splines.
%icp = -[1 3 5 7 8 9]'; % Modes. Verify as design is altered.
%icp = -1;
icp = 0;

% Initialization.
s = createStructure (Nef,Net,Nev,Ner,Nen,Ned,Neh,Neb,Nwel);

s.zdamp = 0; % 0.008; % Modal structural damping ratio.
                      % or
s.adamp = 0;          % Mass-proportional damping factor.
s.bdamp = 0.0002;     % Stiffness-proportional damping factor.

load 'ClCdCm.txt'
Nfoils = (size(ClCdCm,2) - 1)/3;
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
load('TjaereborgAeroProfile.txt');
load('TjaereborgBladeStructure.txt');
AP = TjaereborgAeroProfile;
BS = TjaereborgBladeStructure;

Nain = size(AP,1);      % Number of rows in the input files.
Nsin = size(BS,1);

Lb   = 29.09   /length; % Blade length from pitch bearing flange to tip.

% Get the radial coordinates at which the aero and structural data
% is given.
rbain = Lb*AP(:,1)/length;
rbsin = Lb*BS(:,1)/length;

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
ksei    = ABstruct*BS;
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

cgy(12) = 0;  % Needed for stability, for some reason.

% Aerodynamic properties.
kai     = ABaero*AP;
ai      = Saero*kai;             % Values at elements.
ani     = Snaero*kai;            % Values at nodes.

chord   = ai(:,2)/(1000*length);
xia     = ai(:,3)*pi/180;        % Aero twist, rad.  Positive is LE into wind.
xpc     = ai(:,4)/100;           % Origin of section CS aft of the LE,
                                 % fraction of chord.
tca     = ai(:,5)/100;           % t/c.

% Element coordinates, pitch CS.
xbaero  = [rbel ai(:,8)/(1000*length) ai(:,7)/(1000*length)].';

acent   = ai(:,9)/100;           % Aero center aft of LE, fraction of chord.

% Nodal coordinates, pitch CS.
xbnaero = [rbnod ani(:,8)/(1000*length) ani(:,7)/(1000*length)].'; 
Lel     = norm(xbnaero(:,2:Neb+1) - xbnaero(:,1:Neb),2,'cols').';

% Compute foil weights from t/c.
load('-ascii','tcfoils_Tjaereborg.txt');
foilwt = getFoilWeights (tcfoils_Tjaereborg,tca);

% Put values into the data structures.
for ib = 1:s.Nb

   s.blade(ib).Lb      = Lb;
   s.blade(ib).Lel(:)  = Lel;
   s.blade(ib).xis(:)  = xis;
   s.blade(ib).zeta    = s.zdamp;

   % Store the element twist values in the outboard node.  This is because
   % Node 1 will take on a different meaning as the reference node.
   idof = [10:6:6*(Neb+1)-2].';
   s.blade(ib).Pn_B(idof)    = -xia + xis;
   s.blade(ib).Pn_B(4)       = -xia(1) + xis(1);

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

load('-binary','aoas_Tjaereborg.bin');
load('-binary','kfoils_Tjaereborg.bin');
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
phi   = 0.0    *pi/180; % The root cone angle, for the blade coordinate
                        % system relative to the hub coordinate system.
Ld    = 5.1    /length; % Distance from rear bearing to hub center.
Lbrg  = 2.68   /length; % Bearing-to-bearing length.

Lh    = 1.46   /length; % Hub radius.
mh    = 22100  /mass;   % Hub mass.
kfact = 5;              % Stiffness factor: kfact*k(blade root) = k(hub).

mdrv  = 21000  /mass;   % Driveshaft mass.
Idrv  = 8e5    /(mass*(length^2)); % Effective drv/gear/gen inertia.

rho   = 7850   /(mass/(length^3));
rhoA  = mdrv/Ld;
rhoJ  = 0.05*Idrv/Ld;   % A small nonzero value.  Inertia is accounted for
                        % in the last element, below.
A     = rhoA/rho;
EA    = A*2e11 /(force/(length^2));
Kbend = 1e9    /(force*length);
Ktors = 1.25e8 /(force*length);  % Modified slightly to give the correct
                                 % f_shaft = 0.75 Hz.
EI    = Kbend*Lbrg/3;  % Roark Table 8.1 case 3e.
GJ    = Ktors*Ld;

% The driveshaft.  See modelling guidelines of Snel and Schepers 
% (1995) p 322.
% x/Ld D xis(deg) rhoA rhoJ EA EIyy EIzz GJ
drvprops = [ ...
0.0  1.0/length  0.0  rhoA  rhoJ  EA  EI  EI  GJ; ...
0.5  1.0/length  0.0  rhoA  rhoJ  EA  EI  EI  GJ; ...
1.0  1.0/length  0.0  rhoA  rhoJ  EA  EI  EI  GJ];

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

Dh         = sei(:,2);
xis        = sei(:,3);
rhosA      = sei(:,4);
rhosJ      = sei(:,5);
EA         = sei(:,6);
EIyy       = sei(:,7);
EIzz       = sei(:,8);
GJ         = sei(:,9);

% Additional generator inertia.
igen = 1;
rhosJ(igen) = rhosJ(igen) + Idrv/Lel(igen);

% Fill in driveshaft body entries.
s.driveshaft.Ld             = Ld;
s.driveshaft.Lbrg           = Lbrg;
s.driveshaft.Lh             = Lh;
s.driveshaft.phi            = phi;
s.driveshaft.khfac          = kfact;
s.driveshaft.mh             = mh;

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
Lv    = 3.1    /length; % Vertical distance from yaw bearing to shaft axis.
Lr    = 1.72   /length; % Distance along shaft from yaw center to rear bearing.
Ln    = Lbrg + Lr;      % Distance along shaft from yaw center to front bearing.

delta = 3      *pi/180; % Shaft tilt angle, rad.

mgen  = 80600  /mass;   % Generator mass.
mgear = 16200  /mass;   % Gearbox mass.
mbp   = 21600  /mass;   % Bedplate mass.
mnac  = 34000  /mass;   % Nacelle mass.
maux  = 65400  /mass;   % Auxiliary components' mass.
mbear = 9400   /mass;   % Bearing (x2) mass.

mvr   = mgen + mgear + maux + 0.5*(mbp + mnac);
mnose = mbear + 0.5*(mbp + mnac);

rhoA1 = mvr/(Lv+Lr);
rhoA3 = mnose/Lbrg;
rhoA2 = 0.5*(rhoA1 + rhoA3);

rhoJ  = 3e4/(mass*length);
EA = 2e11/force;                % Make the nacelle elements fairly stiff.
EI = 1e11/(force*(length^2));
GJ = 2e11/(force*(length^2));

% The nacelle turret and nose are parameterized beginning at the 
% yaw bearing, with elements arranged vertically up to the 
% intersection with the driveshaft axis, then horizontally (with
% tilt delta) to the nose.  Nev elements are vertical, Ner to
% the rear bearing (nose flange), and Nen to the front bearing.
% x/Lbrg D xis(deg) rhoA rhoJ EA EIyy EIzz GJ
nacprops = [ ...
0.0         3  0  rhoA1  rhoJ  EA  EI  EI  GJ; ...
Lv/(Lv+Ln)  3  0  rhoA2  rhoJ  EA  EI  EI  GJ; ...
1.0         3  0  rhoA3  rhoJ  EA  EI  EI  GJ];

NevNen     =  Nev + Ner + Nen;
znorm      =  nacprops(:,1);
zin        =  znorm*(Lv + Ln);
zn         =  [([0:Nev].')*Lv/Nev;      ...
               Lv + ([1:Ner].')*Lr/Ner; ...
               Lv + Lr + ([1:Nen].')*(Ln - Lr)/Nen];
ze         =  0.5*(zn(1:NevNen) + zn(2:NevNen+1));
Lel        =  -zn(1:NevNen) + zn(2:NevNen+1);

sei = interp1 (zin,nacprops,ze);

Dh         = sei(:,2);
xis        = sei(:,3);
rhosA      = sei(:,4);
rhosJ      = sei(:,5);
EA         = sei(:,6);
EIyy       = sei(:,7);
EIzz       = sei(:,8);
GJ         = sei(:,9);

s.nacelle.Lv     = Lv;
s.nacelle.Lr     = Lr;
s.nacelle.Ln     = Ln;
s.nacelle.delta  = delta;

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
% The tower here represents the upper half of the actual
% tower.  The foundation represents the lower half.
Lt   = 28       /length;  % Tower height (half).
tt   = 0.25     /length;  % Tower thickness.
rhot = 2600     /(mass/length^3); % Tower density.
Etow = 4.7e10   /(force/(length^2)); % To be calibrated. ft = 0.81 Hz.
Cdt0 = 0.65;

Dt = [7.25;4.75;4.25]/length;
At = 0.25*pi*(Dt.^2 - (Dt - 2*tt).^2);
rhoA = rhot*At;
Jt = (pi/32)*(Dt.^4 - (Dt - 2*tt).^4);
rhoJ = rhot*Jt;
EA = Etow*At;
II = 0.5*Jt;
EI = Etow*II;
GJ = (Etow/2.6)*Jt;

% x/Lt Dt xis(deg) rhoA rhoJ EA EIyy EIzz GJ
towprops = [ ...
0.00  Dt(2)  0  rhoA(2)  rhoJ(2)  EA(2)  EI(2)  EI(2)  GJ(2); ...
0.50  Dt(2)  0  rhoA(2)  rhoJ(2)  EA(2)  EI(2)  EI(2)  GJ(2); ...
1.00  Dt(3)  0  rhoA(3)  rhoJ(3)  EA(3)  EI(3)  EI(3)  GJ(3)];

% Avoid numerically perfect symmetry.
towprops(:,7) = (1 + 1e-6)*towprops(:,7);

% Use linear interpolation.
znorm      =  towprops(:,1);
zin        =  znorm*Lt;
zn         =  ([0:Net].')*Lt/Net;
ze         =  0.5*(zn(1:Net) + zn(2:Net+1));
Lel        = -zn(1:Net) + zn(2:Net+1);

sei = interp1 (zin,towprops,ze);

Dh         = sei(:,2);
xis        = sei(:,3);
rhosA      = sei(:,4);
rhosJ      = sei(:,5);
EA         = sei(:,6);
EIyy       = sei(:,7);
EIzz       = sei(:,8);
GJ         = sei(:,9);

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
   s.tower.Pn_B(idof+3) = zn(inod) - zn(1);
end

%==========================================================
% Foundation.
%==========================================================
Lf = 28      /length;  % Half tower height.

wel = 0;           % Must match length of Nwel at top of input file.

fnprops = [ ...
0.00  Dt(1)  0  rhoA(1)  rhoJ(1)  EA(1)  EI(1)  EI(1)  GJ(1); ...
0.50  Dt(2)  0  rhoA(2)  rhoJ(2)  EA(2)  EI(2)  EI(2)  GJ(2); ...
1.00  Dt(2)  0  rhoA(2)  rhoJ(2)  EA(2)  EI(2)  EI(2)  GJ(2)];

% Avoid numerically perfect symmetry.
fnprops(:,7) = (1 + 1e-6)*fnprops(:,7);

% Use linear interpolation.
znorm      =  fnprops(:,1);
zin        =  znorm*Lf;
zn         =  ([0:Nef].')*Lf/Nef;
ze         =  0.5*(zn(1:Nef) + zn(2:Nef+1));
Lel        = -zn(1:Nef) + zn(2:Nef+1);

sei = interp1 (zin,fnprops,ze);

Dh         = sei(:,2);
xis        = sei(:,3);
rhosA      = sei(:,4);
rhosJ      = sei(:,5);
EA         = sei(:,6);
EIyy       = sei(:,7);
EIzz       = sei(:,8);
GJ         = sei(:,9);

% Soil spring properties.  Adapted from ORT foundation.



%                                  --------------------
scale = 0.0;  % <================= | No soil springs. |
%                                  --------------------



dkxy = scale*[1.e8; 1.e8; 1.e8;  1.e8; 1.e8]/(force/length);

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



