function [M,G,H,K,C] = assembleElementsNL (s,idofs,qq,dqqdt,Pn)
%
% Assemble the nonlinear structural terms.  M multiplies d2q/dt2, G
% multiplies dq/dt, and H and K are vectors.
%
% Version:        Changes:
% --------        -------------
% 03.02.2018      Original code, based on assembleElements.m.
%
% Version:        Verification:
% --------        -------------
% 03.02.2018      
%
% Inputs:
% -------
% s               : Wind turbine data structure.
% idofs           : DOFs references of slave nodes.
% qq...d2qqdt2    : Displacements, velocities, accelerations.
% Pn              : Undeformed nodal positions.
%
% Outputs:
% --------
% M,G,H,K         : Assembled mass and gyroscopic matrix, centrifugal
%                   and stiffness vectors.
% C               : Damping matrix.

%'assembleElementsNL'

Ndof = size(qq,1);  % May include joint angles etc appended to DOFs.
                    % These extra DOFs will remain zero in the matrices.

M = spalloc(Ndof,Ndof,36*Ndof);
G = spalloc(Ndof,Ndof,36*Ndof);
H = zeros(Ndof,1);
K = zeros(Ndof,1);
C = spalloc(Ndof,Ndof,36*Ndof);

for ibod = 1:7

   if (ibod == 1)
      idref = idofs(1);
      Nel   = s.foundation.Nel;
      conns = s.foundation.conn;
      mes   = s.foundation.me_s;
      kes   = s.foundation.ke_s;
   elseif (ibod == 2)
      idref = idofs(2);
      Nel   = s.tower.Nel;
      conns = s.tower.conn;
      mes   = s.tower.me_s;
      kes   = s.tower.ke_s;
   elseif (ibod == 3)
      idref = idofs(3);
      Nel   = s.nacelle.Nel;
      conns = s.nacelle.conn;
      mes   = s.nacelle.me_s;
      kes   = s.nacelle.ke_s;
   elseif (ibod == 4)
      idref = idofs(4);
      Nel   = s.driveshaft.Nel;
      conns = s.driveshaft.conn;
      mes   = s.driveshaft.me_s;
      kes   = s.driveshaft.ke_s;
   elseif (ibod == 5)
      idref = idofs(6);
      Nel   = s.blade(1).Nel;
      conns = s.blade(1).conn;
      mes   = s.blade(1).me_s;
      kes   = s.blade(1).ke_s;
   elseif (ibod == 6)
      idref = idofs(7);
      Nel   = s.blade(2).Nel;
      conns = s.blade(2).conn;
      mes   = s.blade(2).me_s;
      kes   = s.blade(2).ke_s;
   elseif (ibod == 7)
      idref = idofs(8);
      Nel   = s.blade(3).Nel;
      conns = s.blade(3).conn;
      mes   = s.blade(3).me_s;
      kes   = s.blade(3).ke_s;
   end

   for iel = 1:Nel

      ic12 = 12*(iel-1);

      conn  = conns(:,iel);
      rdof  = idref + 6*(conn(1)-1);
      n1dof = idref + 6*(conn(2)-1);
      n2dof = idref + 6*(conn(3)-1);

      qB  = qq(rdof+[1:6]);
      qn1 = qq(n1dof+[1:6]);
      qn2 = qq(n2dof+[1:6]);
      PB  = Pn(rdof+[1:6]);
      Pn1 = Pn(n1dof+[1:6]);
      Pn2 = Pn(n2dof+[1:6]);

      dqdt   = dqqdt([rdof+[1:6] n1dof+[1:6] n2dof+[1:6]]);

      if (n1dof == rdof)

         qn1 = zeros(6,1);
         dqdt(7:12) = zeros(6,1);
         Pn1(1:3) = zeros(3,1);

         % For purposes of computing Qu's, give the reference node the
         % same undeformed orientation as the second node.
         Pn1(4:6) = Pn(n1dof+6+[4:6]);  

      end

      if (n2dof == rdof)
         qn2 = zeros(6,1);
         dqdt(13:18) = zeros(6,1);
         Pn2(1:3) = zeros(3,1);
         Pn2(4:6) = Pn(n2dof+6+[4:6]);
      end

      if (iel > 1)
         if (conns(3,iel-1) == conn(2))
            % Avoid recomputing these nodal matrices where practical.
            Qu1   =   Qu2;
            dQu1  =   dQu2;
            Qu2   =   Qunod (qn2,qB,Pn2,PB);
            dQu2  =   dQudq (qn2,qB,Pn2,PB);
         else
            Qu1   =   Qunod (qn1,qB,Pn1,PB);
            Qu2   =   Qunod (qn2,qB,Pn2,PB);
            dQu1  =   dQudq (qn1,qB,Pn1,PB);
            dQu2  =   dQudq (qn2,qB,Pn2,PB);
         end
      else
         Qu1   =   Qunod (qn1,qB,Pn1,PB);
         Qu2   =   Qunod (qn2,qB,Pn2,PB);
         dQu1  =   dQudq (qn1,qB,Pn1,PB);
         dQu2  =   dQudq (qn2,qB,Pn2,PB);
      end

      [xe,TsB] = elementCSFromNodes (qn1,qn2,Pn1,Pn2);
      dTsB = derivElementCS (qn1,qn2,Pn1,Pn2,TsB);
      mu = getMu (qn1,qn2,Pn1,Pn2,TsB);
      dmu = dmudq (qn1,qn2,Pn1,Pn2,TsB,dTsB);

      [mm,gg,hh,kk] = buildElementNL (mes(:,ic12+[1:12]),kes(:,ic12+[1:12]), ...
                                      dqdt,Qu1,Qu2,dQu1,dQu2,mu,dmu,TsB,dTsB);

      cc = dampingCNL (s.adamp*mes(:,ic12+[1:12])  ...
         +             s.bdamp*kes(:,ic12+[1:12]),dmu);

      M(rdof+[1:6],rdof+[1:6]) = M(rdof+[1:6],rdof+[1:6]) + mm(1:6,1:6);
      G(rdof+[1:6],rdof+[1:6]) = G(rdof+[1:6],rdof+[1:6]) + gg(1:6,1:6);
      H(rdof+[1:6]) = H(rdof+[1:6]) + hh(1:6);
      K(rdof+[1:6]) = K(rdof+[1:6]) + kk(1:6);
      C(rdof+[1:6],rdof+[1:6]) = C(rdof+[1:6],rdof+[1:6]) + cc(1:6,1:6);

      if (n1dof ~= rdof)

         % If node 1 is the reference node, then these nodal DOFs represent the
         % pose of the body in the global coordinate system, which is different
         % than the other nodal DOFs, representing the displacement and
         % transformation matrix T_n^n0 with respect to the undeformed shape.
         % That means that the n1dof entries do not exist, and should not be
         % assigned.  Otherwise, assign them here.

         M( rdof+[1:6],n1dof+[1:6]) = M( rdof+[1:6],n1dof+[1:6]) + mm( 1:6 , 7:12);
         M(n1dof+[1:6], rdof+[1:6]) = M(n1dof+[1:6], rdof+[1:6]) + mm( 7:12, 1:6);
         M(n1dof+[1:6],n1dof+[1:6]) = M(n1dof+[1:6],n1dof+[1:6]) + mm( 7:12, 7:12);

         G( rdof+[1:6],n1dof+[1:6]) = G( rdof+[1:6],n1dof+[1:6]) + gg( 1:6 , 7:12);
         G(n1dof+[1:6], rdof+[1:6]) = G(n1dof+[1:6], rdof+[1:6]) + gg( 7:12, 1:6);
         G(n1dof+[1:6],n1dof+[1:6]) = G(n1dof+[1:6],n1dof+[1:6]) + gg( 7:12, 7:12);

         H(n1dof+[1:6]) = H(n1dof+[1:6]) + hh(7:12);

         K(n1dof+[1:6]) = K(n1dof+[1:6]) + kk(7:12);

         C( rdof+[1:6],n1dof+[1:6]) = C( rdof+[1:6],n1dof+[1:6]) + cc( 1:6 , 7:12);
         C(n1dof+[1:6], rdof+[1:6]) = C(n1dof+[1:6], rdof+[1:6]) + cc( 7:12, 1:6);
         C(n1dof+[1:6],n1dof+[1:6]) = C(n1dof+[1:6],n1dof+[1:6]) + cc( 7:12, 7:12);

      end

      if (n2dof ~= rdof)

         M( rdof+[1:6],n2dof+[1:6]) = M( rdof+[1:6],n2dof+[1:6]) + mm( 1:6 ,13:18);
         M(n2dof+[1:6], rdof+[1:6]) = M(n2dof+[1:6], rdof+[1:6]) + mm(13:18, 1:6);
         M(n2dof+[1:6],n2dof+[1:6]) = M(n2dof+[1:6],n2dof+[1:6]) + mm(13:18,13:18);

         G( rdof+[1:6],n2dof+[1:6]) = G( rdof+[1:6],n2dof+[1:6]) + gg( 1:6 ,13:18);
         G(n2dof+[1:6], rdof+[1:6]) = G(n2dof+[1:6], rdof+[1:6]) + gg(13:18, 1:6);
         G(n2dof+[1:6],n2dof+[1:6]) = G(n2dof+[1:6],n2dof+[1:6]) + gg(13:18,13:18);

         H(n2dof+[1:6]) = H(n2dof+[1:6]) + hh(13:18);

         K(n2dof+[1:6]) = K(n2dof+[1:6]) + kk(13:18);

         C( rdof+[1:6],n2dof+[1:6]) = C( rdof+[1:6],n2dof+[1:6]) + cc( 1:6 ,13:18);
         C(n2dof+[1:6], rdof+[1:6]) = C(n2dof+[1:6], rdof+[1:6]) + cc(13:18, 1:6);
         C(n2dof+[1:6],n2dof+[1:6]) = C(n2dof+[1:6],n2dof+[1:6]) + cc(13:18,13:18);


      end

      if ((n1dof ~= rdof) && (n2dof ~= rdof))

         M(n1dof+[1:6],n2dof+[1:6]) = M(n1dof+[1:6],n2dof+[1:6]) + mm( 7:12,13:18);         
         M(n2dof+[1:6],n1dof+[1:6]) = M(n2dof+[1:6],n1dof+[1:6]) + mm(13:18, 7:12);

         G(n1dof+[1:6],n2dof+[1:6]) = G(n1dof+[1:6],n2dof+[1:6]) + gg( 7:12,13:18);
         G(n2dof+[1:6],n1dof+[1:6]) = G(n2dof+[1:6],n1dof+[1:6]) + gg(13:18, 7:12);

         C(n1dof+[1:6],n2dof+[1:6]) = C(n1dof+[1:6],n2dof+[1:6]) + cc( 7:12,13:18);
         C(n2dof+[1:6],n1dof+[1:6]) = C(n2dof+[1:6],n1dof+[1:6]) + cc(13:18, 7:12);

      end

   end

end

