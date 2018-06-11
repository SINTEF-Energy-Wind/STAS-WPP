function [SpsiR,SpsiI] = ...
             periodicWind (Nef,Net,Nev,Ner,Nen,Ned,Neh,Neb,Nmud,Nwater, ...
                           Lf,Lt,Lv,Lr,Ln,Ld,Lh,Lb,                     ...
                           Ns,Nt,Lel,Diaf,Diat,delta,phi,df,W,z0,Vinf,Vzrp)
% 
% This function computes the velocity spectra of wind shear and tower
% shadow, converted to multi-blade coordinates.
%
% For variable-speed, an idea is to use the estimated probability
% distribution of speed to compute a smooth spectrum, rather than at
% only discrete multiples of a constant P.  This is not yet implemented.
%
% Version:        Changes:
% --------        -------------
% 31.08.2015      Original code.
% 03.08.2017      Adapted for complex step derivatives.
%
% Version:        Verification:
% --------        -------------
% 31.08.2015      
% 03.08.2017      
%
% Inputs:
% -------
% Nef             : Number of foundation elements, to the
%                   base of the tower.
% Net             : Number of tower elements.
% Nev             : Number of vertical elements between the
%                   top of the tower and the driveshaft 
%                   axis.
% Ner             : Number of elements between the corner
%                   of the nacelle turret and the rear
%                   bearing, aligned with the driveshaft
%                   axis.
% Nen             : Number of elements in the nacelle turret
%                   and nose, aligned with the driveshaft
%                   axis.
% Ned             : Number of elements in the driveshaft.
% Neh             : Number of elements between the front
%                   bearing and the hub.
% Neb             : Number of elements per blade.
% Nmud            : Number of elements below the mudline.
% Nwater          : Number of elements between mudline and waterline.
% Lf...Lb         : Dimensions defining the turbine model:
%                   Lf = foundation depth
%                   Lt = Yaw bearing height
%                   Lv = Yaw bearing to shaft axis
%                   Lr = Yaw axis to rear bearing along shaft
%                   Ln = Yaw axis to front bearing
%                   Ld = Rear bearing to hub
%                   Lh = Hub center to pitch bearing
%                   Lb = Pitch bearing to blade tip
% Ns              : Number of Fourier terms for wind shear.
% Nt              : Number of Fourier terms for tower shadow. 
% Lel             : Element length.
% Diaf, Diat      : Diameter of cylinder, of foundation and
%                   tower elements.
% delta           : Tilt of the driveshaft.
% phi             : Blade cone angle.
% df              : Frequency bin width.
% W               : Rotational speed (rad/s).
% z0              : Surface roughness length for wind shear.
% Vinf            : Nominal remote windspeed, for wind shear.
% Vzrp            : Nominal windspeed at the rotorplane.  May include 
%                   induction, conservative if it does not.  Used for
%                   tower shadow.
%
% Outputs:
% --------
% Spsi            : A matrix of velocity auto- and cross-spectra, size
%                   max(Ns,Nt)-by-(3*Neb)^2.

Nel = Nef + Net + Nev + Ner + Nen + Ned + Neh + 3 + 3*Neb;
Lbel = Lel(Nel-3*Neb+1:Nel);

r = zeros(Neb,1);
r(1) = Lh + 0.5*Lbel(1);
for iel = 2:Neb
   r(iel) = r(iel-1) + 0.5*(Lbel(iel) + Lbel(iel-1));
end

H = Lf + Lt - sum(Lel(1:Nmud+Nwater)) + Lv + (Lr+Ld)*sin(delta);

ztel = zeros(Nef+Net+1,1);
zz0 = -sum(Lel(1:Nmud));
ztel(1) = zz0 + 0.5*Lel(1);
for iel = 2:Nef+Net
   ztel(iel) = ztel(iel-1) + 0.5*Lel(iel-1) + 0.5*Lel(iel);
end
zhub = Lf + Lt - sum(Lel(1:Nmud)) + Lv + (Lr+Ld)*sin(delta);
ztel(Nef+Net+1) = zhub;

Dtel = zeros(Nef+Net+1,1);
Dtel(1:Nef) = Diaf;
Dtel(Nef+1:Nef+Net) = Diat;
Dtel(Nef+Net+1) = Diat(Net);

zbel = zhub - r(1:Neb);

xbel = Lr + Ld + Lh*sin(delta) + (r(1:Neb) - Lh)*sin(delta + phi);

Dbel = interp1(ztel,Dtel,zbel);

[SpsiRs,SpsiIs] = windShearMBC (Neb,r,H,z0,df,W,Vinf,Ns);
[SpsiRt,SpsiIt] = towerShadowMBC (Neb,r,xbel,Dbel,df,W,Vzrp,Nt);

SpsiR = zeros(max(Ns,Nt),size(SpsiRs,2));
SpsiI = zeros(max(Ns,Nt),size(SpsiRs,2));
SpsiR(1:Ns,:) = SpsiRs;
SpsiR(1:Nt,:) = SpsiR(1:Nt,:) + SpsiRt;
SpsiI(1:Ns,:) = SpsiIs;
SpsiI(1:Nt,:) = SpsiI(1:Nt,:) + SpsiIt;

SpsiR = single(SpsiR);
SpsiI = single(SpsiI);
