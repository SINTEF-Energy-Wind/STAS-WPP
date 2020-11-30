function [M,dM,dMg,G,dG,dGd,H,dH,dHd,C,dC,K,dK] = ...
               assembleElements (linFlag,s,qq,dqqdt,d2qqdt2,Pn,grav)
%
% Assemble the structural matrices.
%
% Version:        Changes:
% --------        -------------
% 15.01.2020      Original code, a reformulation of 
%                 assembleElementsLin/NL.  Added gravitational
%                 vector.
%
% Version:        Verification:
% --------        -------------
% 15.01.2020      Returns the same results as assembleElementsLin,
%                 discounting gravity.
%
% Inputs:
% -------
% linFlag         : = 1 for linearization.
% s               : Wind turbine data structure.
% idofs           : DOFs references of slave nodes.
% qq,dqqdt,d2qqdt2: Displacements, velocities, accelerations.
% Pn              : Undeformed nodal positions.
% grav            : Gravitational acceleration vector.
%
% Outputs:
% --------
% M,G,H,K         : Assembled mass and gyroscopic matrix, centrifugal
%                   and stiffness vectors.
% C               : Damping matrix.
% dM...dK         : dM/dq*d2q/dt2 etc.

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

Ndof = size(qq,1);  % May include joint angles etc appended to DOFs.
                    % These extra DOFs will remain zero in the matrices.

% Could make these full to trade memory for computation time, but
% the gains seemed negligible.
M   = sparse(Ndof,Ndof);  
G   = sparse(Ndof,Ndof);
H   = zeros(Ndof,1);
K   = zeros(Ndof,1);
C   = sparse(Ndof,Ndof);

dM  = sparse(Ndof,Ndof);
dMg = sparse(Ndof,Ndof);
dG  = sparse(Ndof,Ndof);
dGd = sparse(Ndof,Ndof);
dH  = sparse(Ndof,Ndof);
dHd = sparse(Ndof,Ndof);
dC  = sparse(Ndof,Ndof);
dK  = sparse(Ndof,Ndof);

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

      dqdt = dqqdt([rdof+[1:6] n1dof+[1:6] n2dof+[1:6]]);
      d2qdt2 = d2qqdt2([rdof+[1:6] n1dof+[1:6] n2dof+[1:6]]);
      gvec = [grav;zeros(15,1)];

      if (n1dof == rdof)

         qn1 = zeros(6,1);
         dqdt(7:12) = zeros(6,1);
         d2qdt2(7:12) = zeros(6,1);
         Pn1(1:3) = zeros(3,1);

         % For purposes of computing Qu's, give the reference node the
         % same undeformed orientation as the second node.
         Pn1(4:6) = Pn(n1dof+6+[4:6]);  

      end

      if (n2dof == rdof)
         qn2 = zeros(6,1);
         dqdt(13:18) = zeros(6,1);
         d2qdt2(13:18) = zeros(6,1);
         Pn2(1:3) = zeros(3,1);
         Pn2(4:6) = Pn(n2dof+6+[4:6]);
      end

      if (iel > 1)
         if (conns(3,iel-1) == conn(2))
            % Avoid recomputing these nodal matrices where practical.
            Qu1   = Qu2;
            dQu1  = dQu2;
            d2Qu1 = d2Qu2;
            Qu2   =   Qunod (qn2,qB,Pn2,PB);
            dQu2  =   dQudq (qn2,qB,Pn2,PB);
            if (linFlag == 1)
               d2Qu2 = d2Qudq2 (qn2,qB,Pn2,PB);
            end
         else
            Qu1   =   Qunod (qn1,qB,Pn1,PB);
            Qu2   =   Qunod (qn2,qB,Pn2,PB);
            dQu1  =   dQudq (qn1,qB,Pn1,PB);
            dQu2  =   dQudq (qn2,qB,Pn2,PB);
            if (linFlag == 1)
               d2Qu1 = d2Qudq2 (qn1,qB,Pn1,PB);
               d2Qu2 = d2Qudq2 (qn2,qB,Pn2,PB);
            end
         end
      else
         Qu1   =   Qunod (qn1,qB,Pn1,PB);
         Qu2   =   Qunod (qn2,qB,Pn2,PB);
         dQu1  =   dQudq (qn1,qB,Pn1,PB);
         dQu2  =   dQudq (qn2,qB,Pn2,PB);
         if (linFlag == 1)
            d2Qu1 = d2Qudq2 (qn1,qB,Pn1,PB);
            d2Qu2 = d2Qudq2 (qn2,qB,Pn2,PB);
         end
      end

      if (linFlag ~= 1)
         d2Qu1 = sparse(6,12*12*12);
         d2Qu2 = sparse(6,12*12*12);
      end

      [mu,dmu,d2mu,TsB,dTsB,d2TsB] = ...
                   prepareElementInput (linFlag,qB,qn1,qn2,PB,Pn1,Pn2);

      [mm,dmm,dmmg,gg,dgg,dggd,hv,dhv,dhvd,kv,dkv] =              ...
              buildElement (linFlag,mes(:,ic12+[1:12]),      ...
                            kes(:,ic12+[1:12]),              ...
                            dqdt,d2qdt2,                     ...
                            Qu1,Qu2,dQu1,dQu2,d2Qu1,d2Qu2,   ...
                            mu,dmu,d2mu,TsB,dTsB,d2TsB,gvec);

      [cc,dcc] = dampingC (linFlag,s.adamp*mes(:,ic12+[1:12])  ...
               +           s.bdamp*kes(:,ic12+[1:12]),         ...
                           dqdt,dmu,d2mu);

      M(rdof+[1:6],rdof+[1:6])   = M(rdof+[1:6],rdof+[1:6])   + mm(1:6,1:6);
      G(rdof+[1:6],rdof+[1:6])   = G(rdof+[1:6],rdof+[1:6])   + gg(1:6,1:6);
      H(rdof+[1:6])              = H(rdof+[1:6])              + hv(1:6);
      C(rdof+[1:6],rdof+[1:6])   = C(rdof+[1:6],rdof+[1:6])   + cc(1:6,1:6);
      K(rdof+[1:6])              = K(rdof+[1:6])              + kv(1:6);

      if (linFlag == 1)
         dM(rdof+[1:6],rdof+[1:6])  = dM(rdof+[1:6],rdof+[1:6])  + dmm(1:6,1:6);
         dMg(rdof+[1:6],rdof+[1:6]) = dMg(rdof+[1:6],rdof+[1:6]) + dmmg(1:6,1:6);
         dG(rdof+[1:6],rdof+[1:6])  = dG(rdof+[1:6],rdof+[1:6])  + dgg(1:6,1:6);
         dGd(rdof+[1:6],rdof+[1:6]) = dGd(rdof+[1:6],rdof+[1:6]) + dggd(1:6,1:6);
         dH(rdof+[1:6],rdof+[1:6])  = dH(rdof+[1:6],rdof+[1:6])  + dhv(1:6,1:6);
         dHd(rdof+[1:6],rdof+[1:6]) = dHd(rdof+[1:6],rdof+[1:6]) + dhvd(1:6,1:6);
         dC(rdof+[1:6],rdof+[1:6])  = dC(rdof+[1:6],rdof+[1:6])  + dcc(1:6,1:6);
         dK(rdof+[1:6],rdof+[1:6])  = dK(rdof+[1:6],rdof+[1:6])  + dkv(1:6,1:6);
      end

      if (n1dof ~= rdof)

         % If node 1 is the reference node, then these nodal DOFs represent the
         % pose of the body in the global coordinate system, which is different
         % than the other nodal DOFs, representing the displacement and
         % transformation matrix T_n^n0 with respect to the undeformed shape.
         % That means that the n1dof entries do not exist, and should not be
         % assigned.  Otherwise, assign them here.

         H(n1dof+[1:6]) = H(n1dof+[1:6]) + hv(7:12);
         K(n1dof+[1:6]) = K(n1dof+[1:6]) + kv(7:12);

         M(rdof+[1:6],n1dof+[1:6])   = M(rdof+[1:6],n1dof+[1:6])   + mm(1:6,7:12);
         G(rdof+[1:6],n1dof+[1:6])   = G(rdof+[1:6],n1dof+[1:6])   + gg(1:6,7:12);
         C(rdof+[1:6],n1dof+[1:6])   = C(rdof+[1:6],n1dof+[1:6])   + cc(1:6,7:12);
         M(n1dof+[1:6],rdof+[1:6])   = M(n1dof+[1:6],rdof+[1:6])   + mm(7:12,1:6);
         G(n1dof+[1:6],rdof+[1:6])   = G(n1dof+[1:6],rdof+[1:6])   + gg(7:12,1:6);
         C(n1dof+[1:6],rdof+[1:6])   = C(n1dof+[1:6],rdof+[1:6])   + cc(7:12,1:6);
         M(n1dof+[1:6],n1dof+[1:6])  = M(n1dof+[1:6],n1dof+[1:6])  + mm(7:12,7:12);
         G(n1dof+[1:6],n1dof+[1:6])  = G(n1dof+[1:6],n1dof+[1:6])  + gg(7:12,7:12);
         C(n1dof+[1:6],n1dof+[1:6])  = C(n1dof+[1:6],n1dof+[1:6])  + cc(7:12,7:12);

         if (linFlag == 1)
            dM(rdof+[1:6],n1dof+[1:6])  = dM(rdof+[1:6],n1dof+[1:6])  + dmm(1:6,7:12);
            dMg(rdof+[1:6],n1dof+[1:6]) = dMg(rdof+[1:6],n1dof+[1:6]) + dmmg(1:6,7:12);
            dG(rdof+[1:6],n1dof+[1:6])  = dG(rdof+[1:6],n1dof+[1:6])  + dgg(1:6,7:12);
            dGd(rdof+[1:6],n1dof+[1:6]) = dGd(rdof+[1:6],n1dof+[1:6]) + dggd(1:6,7:12);
            dH(rdof+[1:6],n1dof+[1:6])  = dH(rdof+[1:6],n1dof+[1:6])  + dhv(1:6,7:12);
            dHd(rdof+[1:6],n1dof+[1:6]) = dHd(rdof+[1:6],n1dof+[1:6]) + dhvd(1:6,7:12);
            dC(rdof+[1:6],n1dof+[1:6])  = dC(rdof+[1:6],n1dof+[1:6])  + dcc(1:6,7:12);
            dK(rdof+[1:6],n1dof+[1:6])  = dK(rdof+[1:6],n1dof+[1:6])  + dkv(1:6,7:12);
            dM(n1dof+[1:6],rdof+[1:6])  = dM(n1dof+[1:6],rdof+[1:6])  + dmm(7:12,1:6);
            dMg(n1dof+[1:6],rdof+[1:6]) = dMg(n1dof+[1:6],rdof+[1:6]) + dmmg(7:12,1:6);
            dG(n1dof+[1:6],rdof+[1:6])  = dG(n1dof+[1:6],rdof+[1:6])  + dgg(7:12,1:6);
            dGd(n1dof+[1:6],rdof+[1:6]) = dGd(n1dof+[1:6],rdof+[1:6]) + dggd(7:12,1:6);
            dH(n1dof+[1:6],rdof+[1:6])  = dH(n1dof+[1:6],rdof+[1:6])  + dhv(7:12,1:6);
            dHd(n1dof+[1:6],rdof+[1:6]) = dHd(n1dof+[1:6],rdof+[1:6]) + dhvd(7:12,1:6);
            dC(n1dof+[1:6],rdof+[1:6])  = dC(n1dof+[1:6],rdof+[1:6])  + dcc(7:12,1:6);
            dK(n1dof+[1:6],rdof+[1:6])  = dK(n1dof+[1:6],rdof+[1:6])  + dkv(7:12,1:6);
            dM(n1dof+[1:6],n1dof+[1:6]) = dM(n1dof+[1:6],n1dof+[1:6]) + dmm(7:12,7:12);
            dMg(n1dof+[1:6],n1dof+[1:6])= dMg(n1dof+[1:6],n1dof+[1:6])+ dmmg(7:12,7:12);
            dG(n1dof+[1:6],n1dof+[1:6]) = dG(n1dof+[1:6],n1dof+[1:6]) + dgg(7:12,7:12);
            dGd(n1dof+[1:6],n1dof+[1:6])= dGd(n1dof+[1:6],n1dof+[1:6])+ dggd(7:12,7:12);
            dH(n1dof+[1:6],n1dof+[1:6]) = dH(n1dof+[1:6],n1dof+[1:6]) + dhv(7:12,7:12);
            dHd(n1dof+[1:6],n1dof+[1:6])= dHd(n1dof+[1:6],n1dof+[1:6])+ dhvd(7:12,7:12);
            dC(n1dof+[1:6],n1dof+[1:6]) = dC(n1dof+[1:6],n1dof+[1:6]) + dcc(7:12,7:12);
            dK(n1dof+[1:6],n1dof+[1:6]) = dK(n1dof+[1:6],n1dof+[1:6]) + dkv(7:12,7:12);
         end

      end

      if (n2dof ~= rdof)

         H(n2dof+[1:6]) = H(n2dof+[1:6]) + hv(13:18);
         K(n2dof+[1:6]) = K(n2dof+[1:6]) + kv(13:18);

         M(rdof+[1:6],n2dof+[1:6])   = M(rdof+[1:6],n2dof+[1:6])   + mm(1:6,13:18);
         G(rdof+[1:6],n2dof+[1:6])   = G(rdof+[1:6],n2dof+[1:6])   + gg(1:6,13:18);
         C(rdof+[1:6],n2dof+[1:6])   = C(rdof+[1:6],n2dof+[1:6])   + cc(1:6,13:18);
         M(n2dof+[1:6],rdof+[1:6])   = M(n2dof+[1:6],rdof+[1:6])   + mm(13:18,1:6);
         G(n2dof+[1:6],rdof+[1:6])   = G(n2dof+[1:6],rdof+[1:6])   + gg(13:18,1:6);
         C(n2dof+[1:6],rdof+[1:6])   = C(n2dof+[1:6],rdof+[1:6])   + cc(13:18,1:6);
         M(n2dof+[1:6],n2dof+[1:6])  = M(n2dof+[1:6],n2dof+[1:6])  + mm(13:18,13:18);
         G(n2dof+[1:6],n2dof+[1:6])  = G(n2dof+[1:6],n2dof+[1:6])  + gg(13:18,13:18);
         C(n2dof+[1:6],n2dof+[1:6])  = C(n2dof+[1:6],n2dof+[1:6])  + cc(13:18,13:18);

         if (linFlag == 1)
            dM(rdof+[1:6],n2dof+[1:6])  = dM(rdof+[1:6],n2dof+[1:6])  + dmm(1:6,13:18);
            dMg(rdof+[1:6],n2dof+[1:6]) = dMg(rdof+[1:6],n2dof+[1:6]) + dmmg(1:6,13:18);
            dG(rdof+[1:6],n2dof+[1:6])  = dG(rdof+[1:6],n2dof+[1:6])  + dgg(1:6,13:18);
            dGd(rdof+[1:6],n2dof+[1:6]) = dGd(rdof+[1:6],n2dof+[1:6]) + dggd(1:6,13:18);
            dH(rdof+[1:6],n2dof+[1:6])  = dH(rdof+[1:6],n2dof+[1:6])  + dhv(1:6,13:18);
            dHd(rdof+[1:6],n2dof+[1:6]) = dHd(rdof+[1:6],n2dof+[1:6]) + dhvd(1:6,13:18);
            dC(rdof+[1:6],n2dof+[1:6])  = dC(rdof+[1:6],n2dof+[1:6])  + dcc(1:6,13:18);
            dK(rdof+[1:6],n2dof+[1:6])  = dK(rdof+[1:6],n2dof+[1:6])  + dkv(1:6,13:18);
            dM(n2dof+[1:6],rdof+[1:6])  = dM(n2dof+[1:6],rdof+[1:6])  + dmm(13:18,1:6);
            dMg(n2dof+[1:6],rdof+[1:6]) = dMg(n2dof+[1:6],rdof+[1:6]) + dmmg(13:18,1:6);
            dG(n2dof+[1:6],rdof+[1:6])  = dG(n2dof+[1:6],rdof+[1:6])  + dgg(13:18,1:6);
            dGd(n2dof+[1:6],rdof+[1:6]) = dGd(n2dof+[1:6],rdof+[1:6]) + dggd(13:18,1:6);
            dH(n2dof+[1:6],rdof+[1:6])  = dH(n2dof+[1:6],rdof+[1:6])  + dhv(13:18,1:6);
            dHd(n2dof+[1:6],rdof+[1:6]) = dHd(n2dof+[1:6],rdof+[1:6]) + dhvd(13:18,1:6);
            dC(n2dof+[1:6],rdof+[1:6])  = dC(n2dof+[1:6],rdof+[1:6])  + dcc(13:18,1:6);
            dK(n2dof+[1:6],rdof+[1:6])  = dK(n2dof+[1:6],rdof+[1:6])  + dkv(13:18,1:6);
            dM(n2dof+[1:6],n2dof+[1:6]) = dM(n2dof+[1:6],n2dof+[1:6]) + dmm(13:18,13:18);
            dMg(n2dof+[1:6],n2dof+[1:6])= dMg(n2dof+[1:6],n2dof+[1:6])+ dmmg(13:18,13:18);
            dG(n2dof+[1:6],n2dof+[1:6]) = dG(n2dof+[1:6],n2dof+[1:6]) + dgg(13:18,13:18);
            dGd(n2dof+[1:6],n2dof+[1:6])= dGd(n2dof+[1:6],n2dof+[1:6])+ dggd(13:18,13:18);
            dH(n2dof+[1:6],n2dof+[1:6]) = dH(n2dof+[1:6],n2dof+[1:6]) + dhv(13:18,13:18);
            dHd(n2dof+[1:6],n2dof+[1:6])= dHd(n2dof+[1:6],n2dof+[1:6])+ dhvd(13:18,13:18);
            dC(n2dof+[1:6],n2dof+[1:6]) = dC(n2dof+[1:6],n2dof+[1:6]) + dcc(13:18,13:18);
            dK(n2dof+[1:6],n2dof+[1:6]) = dK(n2dof+[1:6],n2dof+[1:6]) + dkv(13:18,13:18);
         end

      end

      if ((n1dof ~= rdof) && (n2dof ~= rdof))

         M(n1dof+[1:6],n2dof+[1:6])   = M(n1dof+[1:6],n2dof+[1:6])   + mm(7:12,13:18);
         G(n1dof+[1:6],n2dof+[1:6])   = G(n1dof+[1:6],n2dof+[1:6])   + gg(7:12,13:18);
         C(n1dof+[1:6],n2dof+[1:6])   = C(n1dof+[1:6],n2dof+[1:6])   + cc(7:12,13:18);
         M(n2dof+[1:6],n1dof+[1:6])   = M(n2dof+[1:6],n1dof+[1:6])   + mm(13:18,7:12);
         G(n2dof+[1:6],n1dof+[1:6])   = G(n2dof+[1:6],n1dof+[1:6])   + gg(13:18,7:12);
         C(n2dof+[1:6],n1dof+[1:6])   = C(n2dof+[1:6],n1dof+[1:6])   + cc(13:18,7:12);

         if (linFlag == 1)
            dM(n1dof+[1:6],n2dof+[1:6])  = dM(n1dof+[1:6],n2dof+[1:6])  + dmm(7:12,13:18);
            dMg(n1dof+[1:6],n2dof+[1:6]) = dMg(n1dof+[1:6],n2dof+[1:6]) + dmmg(7:12,13:18);
            dG(n1dof+[1:6],n2dof+[1:6])  = dG(n1dof+[1:6],n2dof+[1:6])  + dgg(7:12,13:18);
            dGd(n1dof+[1:6],n2dof+[1:6]) = dGd(n1dof+[1:6],n2dof+[1:6]) + dggd(7:12,13:18);
            dH(n1dof+[1:6],n2dof+[1:6])  = dH(n1dof+[1:6],n2dof+[1:6])  + dhv(7:12,13:18);
            dHd(n1dof+[1:6],n2dof+[1:6]) = dHd(n1dof+[1:6],n2dof+[1:6]) + dhvd(7:12,13:18);
            dC(n1dof+[1:6],n2dof+[1:6])  = dC(n1dof+[1:6],n2dof+[1:6])  + dcc(7:12,13:18);
            dK(n1dof+[1:6],n2dof+[1:6])  = dK(n1dof+[1:6],n2dof+[1:6])  + dkv(7:12,13:18);
            dM(n2dof+[1:6],n1dof+[1:6])  = dM(n2dof+[1:6],n1dof+[1:6])  + dmm(13:18,7:12);
            dMg(n2dof+[1:6],n1dof+[1:6]) = dMg(n2dof+[1:6],n1dof+[1:6]) + dmmg(13:18,7:12);
            dG(n2dof+[1:6],n1dof+[1:6])  = dG(n2dof+[1:6],n1dof+[1:6])  + dgg(13:18,7:12);
            dGd(n2dof+[1:6],n1dof+[1:6]) = dGd(n2dof+[1:6],n1dof+[1:6]) + dggd(13:18,7:12);
            dH(n2dof+[1:6],n1dof+[1:6])  = dH(n2dof+[1:6],n1dof+[1:6])  + dhv(13:18,7:12);
            dHd(n2dof+[1:6],n1dof+[1:6]) = dHd(n2dof+[1:6],n1dof+[1:6]) + dhvd(13:18,7:12);
            dC(n2dof+[1:6],n1dof+[1:6])  = dC(n2dof+[1:6],n1dof+[1:6])  + dcc(13:18,7:12);
            dK(n2dof+[1:6],n1dof+[1:6])  = dK(n2dof+[1:6],n1dof+[1:6])  + dkv(13:18,7:12);
         end

      end

   end

end

M   = sparse(M);  
G   = sparse(G);
C   = sparse(C);

if (linFlag == 1)
   dM  = sparse(dM);
   dMg = sparse(dMg);
   dG  = sparse(dG);
   dGd = sparse(dGd);
   dH  = sparse(dH);
   dHd = sparse(dHd);
   dC  = sparse(dC);
   dK  = sparse(dK);
end
