function [Vavg,vg] = periodicWind (s,a,q,P,psis,ppx,ppy,ppz)
%
% Compute the periodic fluctuations in the wind profile due to
% wind shear and tower shadow.
%
% NOTE: Due to the use of pchip and ppval this does not work
% with complex step derivatives.
%
% Version:        Changes:
% --------        -------------
% 09.12.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 09.12.2019      
%
% Inputs:
% -------
% s,a             : Input data structures from STASTurbine.
% q               : Vector of displacements from a static analysis.
% P               : Vector of undeformed nodal positions.
% ppx,ppy,ppz     : Polynomial fits (via pchip), Vx^g,Vy^g,Vz^g as
%                   a function of global z^g coordinate.
%
% Outputs:
% --------
% Vavg            : A 3*Nel vector of circumferential-average
%                   velocity.
% vg              : An Npsi-by-3*Nel matrix of velocity perturbations
%                   in global coordinates.

Npsi = size(psis,1);
cp = cos (psis);
sp = sin (psis);

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

% Get the P vector specifying undeformed nodal positions and
% orientations.
Pin = assemblePin (s);
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);

% Get some transforms, and call the routine from STAS Aeroelastic
% that computes global blade element positions and projections
% onto the rotorplane.
Try = Tn_y;
Tyy0 = TFromTheta (q(idofs(3)+[4:6]));
Ty0g = TFromTheta (P(idofs(3)+[4:6]));
Trg = Ty0g*Tyy0*Try;
Tyg = Ty0g*Tyy0;
[zr,Area,Dp,r,Lp,xeg,xhg,xyg] = BEMprepProjections (s,a,q,P,Try,Trg);

Tdn = zeros (Npsi,9);
Tdn(:,1) =  cp;
Tdn(:,2) =  sp;
Tdn(:,4) = -sp;
Tdn(:,5) =  cp;
Tdn(:,9) = ones(Npsi,1);

Nel = size(xeg,1)/3;
Nel3 = 3*Nel;
Nel2 = 2*Nel;
oNel = ones(Nel,1);
rg = xeg - repmat(xhg,Nel,1);

i2x = [1:2:Nel2-1].';
i2y = [2:2:Nel2].';

i3x = [1:3:Nel3-2].';
i3y = [2:3:Nel3-1].';
i3z = [3:3:Nel3].';

Tgrs = diagmat (Nel,Trg.');

rr = Tgrs*rg;
rmag = sqrt(rr(i3x).^2 + rr(i3y).^2);
rang = atan2 (rr(i3y),rr(i3x));

r0n = zeros(Nel3,1);
r0n(i3y) = -rmag;
r0n(i3z) = rr(i3z);

Trys = diagmat (Nel,Try);
r0y = Trys*r0n;
r0yz = zeros(Nel3,1);
r0yz(i3z) = r0y(i3z);

% Get the yaw-coordinate offset from the tower centerline to blade
% element location.
xhy = (Tyg.')*xhg;
Oyy = (Tyg.')*P(idofs(3)+[1:3]);
aa = r0y(i3x) + xhy(1) - Oyy(1);

% Each element has an associated analytical domain in which the
% potential flow solution for tower shadow is being computed.  For
% this we need the wind speed orthogonal to the tower centerline,
% at the appropriate height.  This height is the z^g coordinate at
% which the blade passes closest to the tower.
Htsg = xhg(3) + r0yz(i3z);
Vag = zeros(Nel3,1);
Vag(i3x) = ppval (ppx,Htsg);
Vag(i3y) = ppval (ppy,Htsg);
Vag(i3z) = ppval (ppz,Htsg);

% The height can also be used to get the tower diameters.
%znodes = [s.foundation.Pn_B(3:6:6*s.foundation.Nnod-3); ...
%          s.tower.Pn_B(9:6:6*s.tower.Nnod-3)];
znodes = [P(idofs(1)+3)+[0;P([idofs(1)+9:6:idofm(2)+3])]; ...
          P(idofs(2)+3)+P([idofs(2)+9:6:idofm(3)+3])];
Dias = [s.foundation.D(1);s.foundation.D;s.tower.D];
ppd = pchip (znodes,Dias);
Ds = ppval (ppd,Htsg);

% The wind coordinate system, associated with each of the above
% z^g points, has its X^w axis aligned with the wind, and the Z^w
% axis aligned with the tower.  Only the X^w and Y^w components
% of wind are considered for the tower shadow calculation.
Tgys = diagmat (Nel,Tyg.');
Vay = Tgys*Vag;
Vamag = sqrt(Vay(i3x).^2 + Vay(i3y).^2);
thw = atan2 (Vay(i3y),Vay(i3x));
ctw = cos(thw);
stw = sin(thw);

% Compute the sigma vector, which is the position within the 
% analysis plane, but rotated to "alpha" coordinates, which are
% aligned with the yaw coordinate system.
rd = rr;  % Assumes input psi = 0 in defining xeg.

Vavg = zeros(Nel3,1);
vg = zeros(Npsi,Nel3);
for iel = 1:Nel

   ir3 = 3*(iel-1);

   rnp = zeros(Npsi,3);
   rnp(:,1) = Tdn(:,1)*rd(ir3+1) + Tdn(:,4)*rd(ir3+2) + Tdn(:,7)*rd(ir3+3);
   rnp(:,2) = Tdn(:,2)*rd(ir3+1) + Tdn(:,5)*rd(ir3+2) + Tdn(:,8)*rd(ir3+3);
   rnp(:,3) = Tdn(:,3)*rd(ir3+1) + Tdn(:,6)*rd(ir3+2) + Tdn(:,9)*rd(ir3+3);

   sy = zeros(Npsi,3);
   sy(:,1) = Try(1,1)*rnp(:,1) + Try(1,2)*rnp(:,2) + Try(1,3)*rnp(:,3) + xhy(1) - Oyy(1);
   sy(:,2) = Try(2,1)*rnp(:,1) + Try(2,2)*rnp(:,2) + Try(2,3)*rnp(:,3);
   sy(:,3) = Try(3,1)*rnp(:,1) + Try(3,2)*rnp(:,2) + Try(3,3)*rnp(:,3) - r0yz(ir3+3);

   pm1 = sign (sy(:,2));
   pm1(pm1==0) = 1;
   sig = [aa(iel)*ones(Npsi,1), ...
          pm1.*sqrt(sy(:,1).^2 + sy(:,2).^2 + sy(:,3).^2 - aa(iel)^2)];

   Taw = [ctw(iel), -stw(iel); ...
          stw(iel),  ctw(iel)];
   xy = zeros(Npsi,2);
   xy(:,1) = Taw(1,1)*sig(:,1) + Taw(1,2)*sig(:,2);
   xy(:,2) = Taw(2,1)*sig(:,1) + Taw(2,2)*sig(:,2);

   vtsw = zeros(Npsi,2);
   den = (xy(:,1).^2 + xy(:,2).^2).^2;
   vtsw(:,1) = -Vamag(iel)*0.5*(Ds(iel)^2)*(xy(:,1).^2 - xy(:,2).^2)./den;
   vtsw(:,2) = -Vamag(iel)*0.5*(Ds(iel)^2)*xy(:,1).*xy(:,2)./den;

   vtsy = zeros(Npsi,3);
   Twa = Taw.';
   vtsy(:,1) = Twa(1,1)*vtsw(:,1) + Twa(1,2)*vtsw(:,2);
   vtsy(:,2) = Twa(2,1)*vtsw(:,1) + Twa(2,2)*vtsw(:,2);

   vtsg = zeros(Npsi,3);
   vtsg(:,1) = Tyg(1,1)*vtsy(:,1) + Tyg(1,2)*vtsy(:,2);
   vtsg(:,2) = Tyg(2,1)*vtsy(:,1) + Tyg(2,2)*vtsy(:,2);
   vtsg(:,3) = Tyg(3,1)*vtsy(:,1) + Tyg(3,2)*vtsy(:,2);

   % Compute the velocity profiles seen by each element, including wind
   % shear.
   xg = zeros(Npsi,3);
   xg(:,1) = Trg(1,1)*rnp(:,1) + Trg(1,2)*rnp(:,2) + Trg(1,3)*rnp(:,3) + xhg(1);
   xg(:,2) = Trg(2,1)*rnp(:,1) + Trg(2,2)*rnp(:,2) + Trg(2,3)*rnp(:,3) + xhg(2);
   xg(:,3) = Trg(3,1)*rnp(:,1) + Trg(3,2)*rnp(:,2) + Trg(3,3)*rnp(:,3) + xhg(3);

   Vwsg = zeros(Npsi,3);
   Vwsg(:,1) = ppval (ppx,xg(:,3));
   Vwsg(:,2) = ppval (ppy,xg(:,3));
   Vwsg(:,3) = ppval (ppz,xg(:,3));

   VV = Vwsg + vtsg;

   Vavg(ir3+[1:3]) = mean(VV).';
   vg(:,ir3+1) = VV(:,1) - Vavg(ir3+1);
   vg(:,ir3+2) = VV(:,2) - Vavg(ir3+2);
   vg(:,ir3+3) = VV(:,3) - Vavg(ir3+3);

%iel
%Vavg(ir3+[1:3])
%[xg Vwsg vtsg]

end

