function [Tar,Trg,TB0g,TsB,TBB0,dTar,dTsB,dTBB0,wg] = ...
                      BEMprepTransforms (s,a,q,dqdt,P,Tas)
%
% Define transforms that are updated on each call to BEMNL or BEMlin
% for dx/dt.  This needs to be called on each timestep if the blades
% are flexible.
%
% Version:        Changes:
% --------        -------------
% 03.03.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 03.03.2018      
%
% Inputs:
% -------
% aoas,kfoils     : Splined airfoil coefficient tables.
%
% Outputs:
% --------
% Tar             : Transform from airfoil to rotorplane coordinates.
%                   3-by-3*Nel.
% Trg             : Transform from rotorplane to global coordinates.
% TsB             : Transform from section to body coordinates.
%                   3-by-3*Nel.
% TBB0            : Transform from deformed to undeformed body coordinates.
%                   3-by-3*Nel, because we don't assume upfront the
%                   association between element and body.
% dTar            : Derivatives of transform from airfoil to rotorplane.
%                   3-by-3*24*Nel.  qy,qB,qn1,qn2.
% dTsB            : Derivatives of setion-to-body transform.
%                   3-by-3*12*Nel.  qn1,qn2.
% dTBB0           : Derivatives of body reference transform.
%                   3-by-9*Nel.
% wg              : 3*Nel, element structural velocity, global coordinates.

%id = tic();

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

Nel = 0;
for ib = 1:s.Nb
   Nel = Nel + s.blade(ib).Nel;
end

% Coordinate transform matrices.
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);
Try = Tn_y;
Tyy0 = TFromTheta (q(idofs(3)+[4:6]));
Ty0g = TFromTheta (P(idofs(3)+[4:6]));
Trg = Ty0g*Tyy0*Try;
Tgr = Trg.';

% Yaw bearing DOFs, needed as a reference for the rotorplane orientation.
qy  = q(idofs(3)+[1:6]);
Py  = P(idofs(3)+[1:6]);

TB0g    = zeros(3,3*Nel);
Tar     = zeros(3,3*Nel);
dTar    = zeros(3,72*Nel);
TsB     = zeros(3,3*Nel);
dTsB    = zeros(3,36*Nel);
TBB0    = zeros(3,3*Nel);
dTBB0   = zeros(3,9*Nel);
wg      = zeros(3*Nel,1);
for ibod = 5:7

   jb3 = 3*(ibod-1);

   if (ibod == 5)
      idref = idofs(6);
      Neb   = s.blade(1).Nel;
      conns = s.blade(1).conn;
      elref = 0;
   elseif (ibod == 6)
      idref = idofs(7);
      Neb   = s.blade(2).Nel;
      conns = s.blade(2).conn;
      elref = elref + s.blade(1).Nel;
   elseif (ibod == 7)
      idref = idofs(8);
      Neb   = s.blade(3).Nel;
      conns = s.blade(3).conn;
      elref = elref + s.blade(2).Nel;
   end

   qB    = q(idref+[1:6]);
   dqBdt = dqdt(idref+[1:6]);
   PB    = P(idref+[1:6]);

   for iel = 1:Neb
%[ibod iel]
      jel = elref + iel;
%id2 = tic();
      jc72  = 72*(jel-1);
      jc36  = 36*(jel-1);
      jc9   =  9*(jel-1);
      jc3   =  3*(jel-1);

      conn  = conns(:,iel);
      rdof  = idref + 6*(conn(1)-1);
      n1dof = idref + 6*(conn(2)-1);
      n2dof = idref + 6*(conn(3)-1);

      qn1    = q(n1dof+[1:6]);
      dqn1dt = dqdt(n1dof+[1:6]);
      qn2    = q(n2dof+[1:6]);
      dqn2dt = dqdt(n2dof+[1:6]);
      Pn1    = P(n1dof+[1:6]);
      Pn2    = P(n2dof+[1:6]);

      if (n1dof == rdof)

         qn1 = zeros(6,1);
         Pn1(1:3) = zeros(3,1);

         dqn1dt = zeros(6,1);

         % Give the reference node the same undeformed orientation
         % as the second node.
         Pn1(4:6) = P(n1dof+6+[4:6]);  

      end

      TB0g(:,jc3+[1:3]) = TFromTheta (PB(4:6));

      [Tar(:,jc3+[1:3]),dTar(:,jc72+[1:72])] = ...
              airfoilToRotor (qy,qB,qn1,qn2,Py,PB,Pn1,Pn2,Try,Tas(:,jc3+[1:3]));
%'a2r'
%toc(id)
      [xe,TsB(:,jc3+[1:3])] = elementCSFromNodes (qn1,qn2,Pn1,Pn2);
      dTsB(:,jc36+[1:36]) = derivElementCS (qn1,qn2,Pn1,Pn2,TsB(:,jc3+[1:3]));
%'elcs'
%toc(id)
      [TBB0(:,jc3+[1:3]),dTBB0(:,jc9+[1:9])] = dTdth (qB(4:6));

      [vng1,ddy] = globalVelocity (qn1,qB,Pn1,PB,dqn1dt,dqBdt); 
      [vng2,ddy] = globalVelocity (qn2,qB,Pn2,PB,dqn2dt,dqBdt); 
      wg(jc3+[1:3]) = 0.5*(vng1(1:3) + vng2(1:3));

%'wg'
%toc(id)

   end

end

%'end BEMprepTransforms'
%toc(id)
