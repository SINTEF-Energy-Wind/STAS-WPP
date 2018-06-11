function [zr,Area,Dp,r,Lp,xeg,xhg,xyg] = ...
                  BEMprepProjections (s,a,q,dqdt,P,Try,Trg)
%
% Define parameters that are updated on each call to BEMNL or BEMlin
% for dx/dt.
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
% zr              : 2*Nel vector, x,y projections of element positions onto
%                   the rotorplane.
% Area            : Projected area of the annulus.
% Dp              : Diameter of the projected outer node on the rotorplane.
% r               : Projected radial coordinate.
% Lp              : Projected element length, not to be confused with Lel.

%id = tic();

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

Nel = 0;
for ib = 1:s.Nb
   Nel = Nel + s.blade(ib).Nel;
end

% Yaw bearing DOFs, needed as a reference for the rotorplane orientation.
qy  = q(idofs(3)+[1:6]);
Py  = P(idofs(3)+[1:6]);

xyg = globalPosition (qy,Py,zeros(6,1),zeros(6,1));

% Hub center location in global coordinates.  This is used as a reference
% for the BEM rotorplane.
xhg = globalPosition (q(idofs(4)+[1:6]),P(idofs(4)+[1:6]), ...
                      q(idofm(6)-6+[1:6]),P(idofm(6)-6+[1:6]));

% Estimate the outer diameter for the induced velocity calculation.
% Take the maximum tip position for each blade, since the analysis
% will fail if the outermost radial coordinate falls outside this.
rtip = zeros(3,1);
for ib = 1:s.Nb
   ic3 = 3*(ib-1);
   xtg = globalPosition (q(idofs(5+ib)+[1:6]),P(idofs(5+ib)+[1:6]), ...
                         q(idofs(5+ib)+6*a.Neb+[1:6]),P(idofs(5+ib)+6*a.Neb+[1:6]));
   xtr = (Trg.')*(xtg - xhg);
   rtip(ib) = sqrt(xtr(1)^2 + xtr(2)^2);
end

%Dp = 2*maxc(rtip(1),maxc(rtip(2),rtip(3)));
Dp = 2*rtip;

zr      = zeros(2*Nel,1);
r       = zeros(Nel,1);
Lp      = zeros(Nel,1);
Area    = zeros(Nel,1);
xeg     = zeros(3*Nel,1);
for ibod = 5:7

   jb3  = 3*(ibod-1);

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
      jc72  = 72*(jel-1);
      jc36  = 36*(jel-1);
      jc9   =  9*(jel-1);
      jc3   =  3*(jel-1);
      jc2   =  2*(jel-1);

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

         % Give the reference node the same undeformed orientation
         % as the second node.
         Pn1(4:6) = P(n1dof+6+[4:6]);  

      end

      [xeg(jc3+[1:3]),xn1,xn2,xer,r(jel),Lp(jel),ddy] = ...
              projectElement (qy,qB,qn1,qn2,Py,PB,Pn1,Pn2,Try,xhg);

      zr(jc2+[1:2]) = xer(1:2);

      Area(jel) = 2*pi*r(jel)*Lp(jel)/s.Nb;
%'end loop'
%toc(id)
   end

end

%'end BEMprepProjections'
%toc(id)