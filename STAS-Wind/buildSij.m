function Sij = buildSij (QSflag,MBCflag,s,a,nnf,df,Vinf,TI,Lu,q,dqdt,P)
%
% Compute the correlation functions or spectra associated with
% turbulence, in either rotating or MBC coordinates.
%
% Version:        Changes:
% --------        -------------
% 12.12.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 12.12.2019      This vectorized version has been compared against
%                 a slower but simpler version with nested "for" 
%                 loops.  Agreement was precise.
%
% Inputs:
% -------
% s,a             : Input data structures from STASTurbine.
% q               : Vector of displacements from a static analysis.
% P               : Vector of undeformed nodal positions.
%
% Outputs:
% --------
% 

Nf = 4*nnf;

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

dt = 1/(Nf*df);
taus = dt*[[0:(Nf/2)-1] [-Nf/2:-1]].';

Vmag = sqrt(Vinf(1)^2 + Vinf(2)^2 + Vinf(3)^2);
sig = TI*Vmag;
sig2 = sig^2;
third = 1/3;
twosiggam = 2*(sig^2)/gamma(third);
oneLu = 1.34*Lu;

iW = idofs(4) + 6;
psis = dqdt(iW)*taus;

cp  = cos (psis);
sp  = sin (psis);

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

Tdn = zeros (Nf,9);
Tdn(:,1) =  cp;
Tdn(:,2) =  sp;
Tdn(:,4) = -sp;
Tdn(:,5) =  cp;
Tdn(:,9) = ones(Nf,1);

Tdn0 = reshape (Tdn(1,:),3,3);

TBpsi = zeros (Nf,9);
if (MBCflag == 1)
   cp2 = cos (psis + 2*pi/3);
   sp2 = sin (psis + 2*pi/3);
   cp3 = cos (psis + 4*pi/3);
   sp3 = sin (psis + 4*pi/3);
   TBpsi(:,1) = third;
   TBpsi(:,2) = 2*third*cp;
   TBpsi(:,3) = 2*third*sp;
   TBpsi(:,4) = third;
   TBpsi(:,5) = 2*third*cp2;
   TBpsi(:,6) = 2*third*sp2;
   TBpsi(:,7) = third;
   TBpsi(:,8) = 2*third*cp3;
   TBpsi(:,9) = 2*third*sp3;
else
   TBpsi(:,1) = 1;
   TBpsi(:,5) = 1;
   TBpsi(:,9) = 1;
end

Nel = size(xeg,1)/3;
Nel3 = 3*Nel;
Nel2 = 2*Nel;
oNel = ones(Nel,1);
Neb = Nel/3;

i3x = [1:3:Nel3-2].';
i3y = [2:3:Nel3-1].';
i3z = [3:3:Nel3].';

ibl = [ones(Neb,1);2*ones(Neb,1);3*ones(Neb,1)];

% Find the constant element coordinates in the driveshaft
% or hub coordinate system.
Tgrs = diagmat (Nel,Trg.');
rg0 = xeg - repmat(xhg,Nel,1);
rd = Tgrs*rg0;

Qij = zeros (Nf,(3*Nel)^2);
for iel1 = 1:Nel

%printf('%5d\n',iel1);
%fflush(stdout);

   ic1 = 9*Nel*(iel1-1);
   i1r3 = 3*(iel1-1);

   ibel1 = mod(iel1-1,Neb)+1;
   itrip1 = [ibel1, ibel1+Neb, ibel1+2*Neb].';
   ic1s = 9*Nel*(itrip1 - 1);

   rnp1 = Tdn0*rd(i1r3+[1:3]);

   Tmbc0 = TBpsi(1,3*(ibl(iel1)-1)+[1:3]);

   for iel2 = 1:Nel

      icol = ic1 + 9*(iel2-1);
      ir3 = 3*(iel2-1);
      ibel2 = mod(iel2-1,Neb)+1;
      itrip2 = [ibel2, ibel2+Neb, ibel2+2*Neb].';

      icols = [ic1s(1) + 9*(itrip2 - 1); ...
               ic1s(2) + 9*(itrip2 - 1); ...
               ic1s(3) + 9*(itrip2 - 1)];

      TmbcT = TBpsi(:,3*(ibl(iel2)-1)+[1:3]);

      r12 = zeros(Nf,3);
      r12(:,1) = Tdn(:,1)*rd(ir3+1) + Tdn(:,4)*rd(ir3+2) + Tdn(:,7)*rd(ir3+3) - rnp1(1);
      r12(:,2) = Tdn(:,2)*rd(ir3+1) + Tdn(:,5)*rd(ir3+2) + Tdn(:,8)*rd(ir3+3) - rnp1(2);
      r12(:,3) = Tdn(:,3)*rd(ir3+1) + Tdn(:,6)*rd(ir3+2) + Tdn(:,9)*rd(ir3+3) - rnp1(3);

      sg = zeros(Nf,3);
      sg(:,1) = Trg(1,1)*r12(:,1) + Trg(1,2)*r12(:,2) + Trg(1,3)*r12(:,3) - Vinf(1)*taus;
      sg(:,2) = Trg(2,1)*r12(:,1) + Trg(2,2)*r12(:,2) + Trg(2,3)*r12(:,3) - Vinf(2)*taus;
      sg(:,3) = Trg(3,1)*r12(:,1) + Trg(3,2)*r12(:,2) + Trg(3,3)*r12(:,3) - Vinf(3)*taus;

      s = sqrt(sg(:,1).^2 + sg(:,2).^2 + sg(:,3).^2);
      sLu = s/oneLu;

      ind = (s > eps^0.25);
      sjs = zeros (Nf,3);
      sjs(ind,1) = sg(ind,1)./s(ind);
      sjs(ind,2) = sg(ind,2)./s(ind);
      sjs(ind,3) = sg(ind,3)./s(ind);
      sjs(~ind,:) = 1;

      Qss = zeros(Nf,1);
      Qss(ind) = twosiggam*((0.5*sLu(ind)).^third).*besselk(third,sLu(ind));
      Qss(~ind) = sig2;

      % This limit as s->0 is actually infinity, but this limit times any of 
      % the s components is zero, so it's safe to set it to zero.
      dQss = zeros(Nf,1);
      dQss(ind) = -(twosiggam/oneLu)*((0.5*sLu(ind)).^third) ...
               .* besselk(-2*third,sLu(ind));
      dQss(~ind) = 0;

      QsdQss2 = Qss + 0.5*s.*dQss;

      QQ = zeros(Nf,9);
      QQ(:,1) = QsdQss2 - 0.5*sjs(:,1).*sg(:,1).*dQss;
      QQ(:,2) = -0.5*sjs(:,1).*sg(:,2).*dQss;
      QQ(:,3) = -0.5*sjs(:,1).*sg(:,3).*dQss;
      QQ(:,4) = QQ(:,2);
      QQ(:,5) = QsdQss2 - 0.5*sjs(:,2).*sg(:,2).*dQss;
      QQ(:,6) = -0.5*sjs(:,2).*sg(:,3).*dQss;
      QQ(:,7) = QQ(:,3);
      QQ(:,8) = QQ(:,6);
      QQ(:,9) = QsdQss2 - 0.5*sjs(:,3).*sg(:,3).*dQss;

      Qij(:,icols(1)+[1:9]) = Qij(:,icols(1)+[1:9]) ...
                             + Tmbc0(1)*QQ.*TmbcT(:,1);
      Qij(:,icols(2)+[1:9]) = Qij(:,icols(2)+[1:9]) ...
                             + Tmbc0(1)*QQ.*TmbcT(:,2);
      Qij(:,icols(3)+[1:9]) = Qij(:,icols(3)+[1:9]) ...
                             + Tmbc0(1)*QQ.*TmbcT(:,3);
      Qij(:,icols(4)+[1:9]) = Qij(:,icols(4)+[1:9]) ...
                             + Tmbc0(2)*QQ.*TmbcT(:,1);
      Qij(:,icols(5)+[1:9]) = Qij(:,icols(5)+[1:9]) ...
                             + Tmbc0(2)*QQ.*TmbcT(:,2);
      Qij(:,icols(6)+[1:9]) = Qij(:,icols(6)+[1:9]) ...
                             + Tmbc0(2)*QQ.*TmbcT(:,3);
      Qij(:,icols(7)+[1:9]) = Qij(:,icols(7)+[1:9]) ...
                             + Tmbc0(3)*QQ.*TmbcT(:,1);
      Qij(:,icols(8)+[1:9]) = Qij(:,icols(8)+[1:9]) ...
                             + Tmbc0(3)*QQ.*TmbcT(:,2);
      Qij(:,icols(9)+[1:9]) = Qij(:,icols(9)+[1:9]) ...
                             + Tmbc0(3)*QQ.*TmbcT(:,3);

   end

end

if (QSflag ~= 1)

   % Perform the Fourier transform to obtain the spectra.  Do this
   % in chunks to save memory.
   M = Nel;
   Nch = size(Qij,2)/M;
   for ich = 1:Nch
      icol = M*(ich-1);
      Qij(:,icol+[1:M]) = dt*fft(Qij(:,icol+[1:M]));
   end

end

Sij = Qij;
