function a = createAero (Nfoils,Neb)
%
% Creates a data structure holding the inputs needed for the
% aerodynamic loads analysis.
%

   Nb       = 3;                      % Hard-coded number of blades.

   a.Nb     = Nb;
   a.Neb    = Neb;                    % Number of aero elements per blade.

   %----------------- Element section properties -----------------
   a.xia    = zeros(Nb*Neb,1);        % Aerodynamic twist, rad.
   a.chord  = zeros(Nb*Neb,1);        % Chord length.
   a.xpc    = zeros(Nb*Neb,1);        % Position of the pitch axis aft of the
                                      % leading edge, in fractions of the chord
                                      % length.
   a.acent  = zeros(Nb*Neb,1);        %   ...           aero center    ...
   a.Lel    = zeros(Nb*Neb,1);        % Length of each aero element, for aero
                                      % load calculation.  Caution, use the
                                      % projected length for momentum balance.
   a.foilwt = sparse(Nfoils,Nb*Neb);  % Weight of each airfoil type.
   a.aoas   = [];
   a.kfoils = [];

   % The zero-lift aoa is used for dynamic stall.  This has to be
   % computed after getting the airfoil coefficient weights.  For
   % "weird" inboard sections that don't quite behave like airfoils,
   % this should be where the projection of the primary lift curve
   % slope intersects the zero-lift axis.
   a.aoaz   = zeros(Nb*Neb,1); 

   a.Pa_B   = zeros(6*Nb*Neb,1);      % Undeformed coordinates of aero elements,
                                      % pitch CS.
   a.Ta_s   = zeros(3,3*Nb*Neb);      % Transform from aero to (aero) section
                                      % coordinates.


   %----------------------- Air properties -----------------------
   a.dens   = 0;                      % Air density.
   a.visc   = 0;                      % Air viscosity.
   a.z0     = 0;                      % Surface roughness length for wind shear.
   
   % Set to positive integers if icp contains the control points
   % for a spline representation of the aerodynamic states. Set
   % to negative integers if icp contains the blade body modes
   % for a modal-slave representation of the aerodynamic states.
   % In the latter case, the selection of modes has to be verified
   % as the turbine design is altered.
   a.icp    = [];

