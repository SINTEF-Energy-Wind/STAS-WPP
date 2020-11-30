function Psia = getPsia (modalAero,s,a,azi,q,qd,P,shape0)

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
Ndj = Ndof + 6;
slv = slaveDOFs (idofs);
vec = [1:Ndj].';
[jnk,ret,jnk2] = partitionMatrix (vec,slv,[]);

cp0 = cos(azi);
sp0 = sin(azi);

[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);
[Tas,ch,Lel,foilwt,aoaz,aoast,xas,yas,iq] = BEMsetup (s,a);

Td_n = [cp0 -sp0 0;sp0 cp0 0;0 0 1];
Try = Tn_y;
Tyy0 = TFromTheta (q(idofs(3)+[4:6]));
Ty0g = TFromTheta (P(idofs(3)+[4:6]));
Trg = Ty0g*Tyy0*Try;

[Tar,Trg,TB0g,TsB,TBB0,dTar,dTsB,dTBB0,wga] = ...
                      BEMprepTransforms (s,a,q,qd,P,Tas);

[zr,Area,Dp,rp,Lp,xeg,xhg,xyg] = ...
         BEMprepProjections (s,a,q,P,Try,Trg);

bsh = bladeModeShape (s,ret,shape0);

% Conversion matrix from aero spline to mode.
if (modalAero == 1)

   abm = a;
   abm.icp = [-1, -3].';  % First two flap modes.
   Psia = aeroPsi (abm,rp,bsh);

else

   Psia = aeroPsi (a,rp,bsh);

end