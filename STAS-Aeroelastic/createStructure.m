function s = createStructure (Nef,Net,Nev,Ner,Nen,Ned,Neh,Neb,Nwel)

   s.Nef        = Nef;
   s.Net        = Net;
   s.Nev        = Nev;
   s.Ner        = Ner;
   s.Nen        = Nen;
   s.Ned        = Ned;
   s.Neh        = Neh;
   s.Neb        = Neb;

   Nb           = 3;
   s.Nb         = Nb;
   s.foundation = createBody(Nef,Nef+1);
   s.tower      = createBody(Net,Net+1);
   s.nacelle    = createBody(Nev+Ner+Nen,Nev+Ner+Nen+1);
   s.driveshaft = createBody(Ned+Neh+3,Ned+Neh+4);
   s.blade(1)   = createBody(Neb,Neb+1);

   % Additional structural parameters, beyond the basic body definition.

   s.foundation.Lf    = 0;       % Height of foundation pile.
   s.foundation.D     = zeros(Nef,1); % Foundation outer diameter at elements.
                                 % Used in the computation of stress.
   s.foundation.rhow  = 0;       % Water density for added mass in mass matrix.
   s.foundation.wel   = zeros(Nwel,1); % List of elements subject to wave loads.
   s.foundation.Nwater= 0;       % Number of submerged elements.
   s.foundation.CmA   = zeros(Nwel,1); % Inertia coefficient times area,
                                 % for added mass in mass matrix.
   s.foundation.Cd    = zeros(Nwel,1); % Drag coefficient.
   s.foundation.k     = zeros(8,Nef+1); % Linearized p-y curve springs: 
                                 % kxg,kyg,kzg,kthzg,cxg,cyg,czg,cthzg
   s.tower.Lt         = 0;       % Height from transition piece to the top
                                 % of the tower.
   s.tower.D          = zeros(Net,1); % Outer tower diameter at elements, for aero
                                 % loads and stress.
   s.tower.Cd         = zeros(Net,1); % Aerodynamic drag coefficient.
   s.nacelle.Lv       = 0;       % Vertical distance from yaw bearing to
                                 % shaft axis.
   s.nacelle.Lr       = 0;       % Distance along shaft from yaw center to
                                 % rear bearing.
   s.nacelle.Ln       = 0;       % Distance along shaft from yaw center to
                                 % front bearing. 
   s.nacelle.delta    = 0;       % Shaft tilt angle, rad.
   s.nacelle.D        = zeros(Nev+Ner+Nen,1); % Nacelle outer diameter at elements.
   s.nacelle.rstat    = 0;       % Effective r of generator stator for Jgen.
   s.nacelle.mgstat   = 0;       % Generator stator mass.
   s.nacelle.Jgstat   = 0;       % Generator stator inertia.
   s.driveshaft.Ld    = 0;       % Distance from rear bearing to hub center.
   s.driveshaft.Lbrg  = 0;       % Bearing-to-bearing length.
   s.driveshaft.Lh    = 0;       % Hub radius.
   s.driveshaft.D     = zeros(Ned+Neh,1); % Driveshaft outer diameter at elements.
   s.driveshaft.phi   = 0;       % Root cone angle of the pitch axis.
   s.driveshaft.mh    = 0;       % Hub mass.
   s.driveshaft.khfac = 0;       % Stiffness factor: khfac*k(blade root) = k(hub).
   s.driveshaft.rrot  = 0;       % Effective r of generator rotor for Jgen.
   s.driveshaft.mgrot = 0;       % Generator rotor mass.
   s.driveshaft.Jgrot = 0;       % Generator rotor inertia.
   s.blade(1).Lb      = 0;       % Blade length from pitch flange to tip.

   for ib = 2:Nb
      s.blade(ib) = s.blade(1);
   end

