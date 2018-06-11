function [F,dFdq,dFdqd] = soilLin (kx,ky,kz,kthz,cx,cy,cz,cthz, ...
                                   PF,qF,dqFdt)
%
% Compute soil stiffness and damping matrices at an operating point.
% We have a soil force F(q,qdot).  This is applied to the equations
% of motion as Qhat(q)*F.  Thus terms in the linear equations of
% motion are
% Ks = -Qhat*(dF/dq), Kf = -(dQhat/dq)*F, Cs = -Qhat*(dF/dqdot).  
% This function returns the operating-point force F, and its
% derivatives.  The Qhat operations are implemented externally.
%
% Version:        Changes:
% --------        -------------
% 13.05.2017      Original code, reformulation of soil.m.
%
% Version:        Verification:
% --------        -------------
% 13.05.2017      
%
% Inputs:
% -------
% kx...cthz       : x,y,z displacement, z rotational spring
%                   stiffness and damping, in global coordinates.
% PF              : Nodal offsets.  The first six are the global
%                   position and orientation of the foundation
%                   reference node.
% qF,dqFdt        : Nodal displacements and velocities.  The first
%                   six are the displacements and orientations of
%                   the foundation reference node.
%
% Outputs:
% --------
% F               : Nodal force vector containing soil stiffness
%                   and damping forces.
% dFdq,dFdqd      : Derivatives wrt q and dq/dt.

Nnod = size(kx,1);
Nd = 6*Nnod;

F     = zeros(Nd,1);
dFdq  = spalloc (Nd,Nd,6*Nd);
dFdqd = spalloc (Nd,Nd,6*Nd);

TF0g = TFromTheta (PF(4:6));
TgF0 = TF0g.';
[TFF0,dTFF0] = dTdth (qF(4:6));
d2TFF0 = d2Tdth2 (qF(4:6),TFF0,dTFF0);
TF0F = TFF0.';

for inod = 1:Nnod

   ic6 = 6*(inod-1);

   Kd = diag([kx(inod) ky(inod) kz(inod)].');
   Cd = diag([cx(inod) cy(inod) cz(inod)].');
   TK = TF0F*Kd;
   TC = TF0F*Cd;

   PB  = PF(1:6);
   qB  = qF(1:6);
   dqB = dqFdt(1:6);

   if (inod == 1)
      % Reference node.
      Pn  = zeros(6,1);
      qn  = zeros(6,1);
      dqn = zeros(6,1);
      Pn(4:6) = PF(6+[4:6]);
   else
      Pn  = PF(ic6+[1:6]);
      qn  = qF(ic6+[1:6]);
      dqn = dqFdt(ic6+[1:6]);
   end

   % Forces due to linear displacements.
   F(ic6+[1:3]) = -TK*(TgF0*qB(1:3)                        ...
                +      TFF0*(Pn(1:3) + qn(1:3)) - Pn(1:3)) ...
                -  TC*(TgF0*dqB(1:3) + TFF0*dqn(1:3));

   for jj = 1:3
      jc3 = 3*(jj-1);
      F(ic6+[1:3]) = F(ic6+[1:3])          ...
                   - TC*dTFF0(:,jc3+[1:3])*(Pn(1:3) + qn(1:3))*dqB(3+jj);
   end

   % Moments due to torsion.
   Tn0F = TFromTheta (Pn(4:6));
   [Tnn0,dTnn0] = dTdth (qn(4:6));
   d2Tnn0 = d2Tdth2 (qn(4:6),Tnn0,dTnn0);
   mat = TFF0*Tn0F*Tnn0;
   x = mat(1,2);
   y = mat(2,2);
   thz = atan2c(y,x);
   F(ic6+[4:6]) = -TF0F*[0 0 kthz(inod)*thz].';

   dT = zeros(3,3);
   for jj = 1:3
      jc3 = 3*(jj-1);
      dT = dT + dTFF0(:,jc3+[1:3])*Tn0F*Tnn0*dqB(3+jj) ...
         +      TFF0*Tn0F*dTnn0(:,jc3+[1:3])*dqn(3+jj);
   end
   dxdt   =  dT(1,2);
   dydt   =  dT(2,2);
   x2y2   =  x^2 + y^2;
   dthdx  = -y/x2y2;
   dthdy  =  x/x2y2;
   dthzdt =  dthdx*dxdt + dthdy*dydt;
   F(ic6+[4:6]) = F(ic6+[4:6]) - TF0F*[0 0 cthz(inod)*dthzdt].';

%del = sqrt(eps);
%yc = y + i*del;
%tc = atan2c (yc,x);
%dtcdy = imag(tc)/del;
%'atan'
%dthdy-dtcdy

   dFdq(ic6+[1:3],1:3) = dFdq(ic6+[1:3],1:3) - TK*TgF0;

   if (inod ~= 1)
      dFdq(ic6+[1:3],ic6+[1:3]) = dFdq(ic6+[1:3],ic6+[1:3]) - TK*TFF0;
   end

   for jj = 1:3
      jc3 = 3*(jj-1);
      jc9 = 9*(jj-1);
      dFdq(ic6+[1:3],3+jj) = dFdq(ic6+[1:3],3+jj)        ...
             - TK*dTFF0(:,jc3+[1:3])*(Pn(1:3) + qn(1:3)) ...
             - TC*dTFF0(:,jc3+[1:3])*dqn(1:3)            ...
             + (dTFF0(:,jc3+[1:3]).')*TFF0*F(ic6+[1:3]);
      if (inod ~= 1)
         dFdq(ic6+[1:3],ic6+[1:3]) = dFdq(ic6+[1:3],ic6+[1:3]) ...
                                   - TC*dTFF0(:,jc3+[1:3])*dqB(3+jj);
      end
      for ii = 1:3
         i3 = 3*(ii-1);
         dFdq(ic6+[1:3],3+jj) = dFdq(ic6+[1:3],3+jj) ...
              - TC*d2TFF0(:,jc9+i3+[1:3])*(Pn(1:3) + qn(1:3))*dqB(3+ii);
      end
   end

   dFdqd(ic6+[1:3],1:3) = dFdqd(ic6+[1:3],1:3) - TC*TgF0;
   if (inod ~= 1)
      dFdqd(ic6+[1:3],ic6+[1:3]) = dFdqd(ic6+[1:3],ic6+[1:3]) - TC*TFF0;
   end
   for jj = 1:3
      jc3 = 3*(jj-1);
      dFdqd(ic6+[1:3],3+jj) = dFdqd(ic6+[1:3],3+jj) ...
                            - TC*dTFF0(:,jc3+[1:3])*(Pn(1:3) + qn(1:3));
   end

   delT = zeros (3,3*6);    % d/d(thB,thn)
   for jj = 1:3
      jc3 = 3*(jj-1);
      delT(:,jc3+[1:3]) = delT(:,jc3+[1:3]) ...
                        + dTFF0(:,jc3+[1:3])*Tn0F*Tnn0;
      delT(:,jc3+9+[1:3]) = delT(:,jc3+9+[1:3]) ...
                          + TFF0*Tn0F*dTnn0(:,jc3+[1:3]);
   end

   ddelT = zeros (3,3*12);  % d/d(thB,thn,dthB,dthn)

   for jj = 1:3
      jc9 = 9*(jj-1);
      jc3 = 3*(jj-1);
      ddelT(:,18+jc3+[1:3]) = ddelT(:,jc3+[1:3]) ...
                            + dTFF0(:,jc3+[1:3])*Tn0F*Tnn0;
      ddelT(:,27+jc3+[1:3]) = ddelT(:,27+jc3+[1:3]) ...
                            + TFF0*Tn0F*dTnn0(:,jc3+[1:3]);
      for ii = 1:3
         i3 = 3*(ii-1);
         ddelT(:,jc3+[1:3]) = ddelT(:,jc3+[1:3])             ...
                + d2TFF0(:,jc9+i3+[1:3])*Tn0F*Tnn0*dqB(3+ii) ...
                + dTFF0(:,jc3+[1:3])*Tn0F*dTnn0(:,i3+[1:3])*dqn(3+ii);
         ddelT(:,9+jc3+[1:3]) = ddelT(:,9+jc3+[1:3])                  ...
                + dTFF0(:,i3+[1:3])*Tn0F*dTnn0(:,jc3+[1:3])*dqB(3+ii) ...
                + TFF0*Tn0F*d2Tnn0(:,jc9+i3+[1:3])*dqn(3+ii);

      end

   end

   dMdthz   = -kthz(inod);
   dMdthzd  = -cthz(inod);
   x2y2     =  x^2 + y^2;
   dthdx    = -y/x2y2;
   dthdy    =  x/x2y2;
   d2thdx2  =  2*x*y/(x2y2^2);
   d2thdy2  = -d2thdx2;
   d2thdxdy = (y^2 - x^2)/(x2y2^2);

   dxdthB = [delT(1,2)  delT(1,5)  delT(1,8)];
   dydthB = [delT(2,2)  delT(2,5)  delT(2,8)];
   dxdthn = [delT(1,11) delT(1,14) delT(1,17)];
   dydthn = [delT(2,11) delT(2,14) delT(2,17)];

   dxddthB  = [ddelT(1,2)  ddelT(1,5)  ddelT(1,8)];
   dyddthB  = [ddelT(2,2)  ddelT(2,5)  ddelT(2,8)];
   dxddthn  = [ddelT(1,11) ddelT(1,14) ddelT(1,17)];
   dyddthn  = [ddelT(2,11) ddelT(2,14) ddelT(2,17)];
   dxddthBd = [ddelT(1,20) ddelT(1,23) ddelT(1,26)];
   dyddthBd = [ddelT(2,20) ddelT(2,23) ddelT(2,26)];
   dxddthnd = [ddelT(1,29) ddelT(1,32) ddelT(1,35)];
   dyddthnd = [ddelT(2,29) ddelT(2,32) ddelT(2,35)];

   dthddx  = d2thdx2*dxdt + d2thdxdy*dydt;
   dthddy  = d2thdxdy*dxdt + d2thdy2*dydt;
   dthddxd = dthdx;
   dthddyd = dthdy;

   dFdq(ic6+[4:6],4:6) = dFdq(ic6+[4:6],4:6)              ...
          + TF0F*[zeros(1,3);zeros(1,3);                  ...
                  dMdthz*(dthdx*dxdthB + dthdy*dydthB)    ...
          +       dMdthzd*(dthddx*dxdthB + dthddy*dydthB) ...
          +       dMdthzd*(dthddxd*dxddthB + dthddyd*dyddthB)];
   for jj = 1:3
      jc3 = 3*(jj-1);
      dFdq(ic6+[4:6],3+jj) = dFdq(ic6+[4:6],3+jj)         ...
          + (dTFF0(:,jc3+[1:3]).')*TFF0*F(ic6+[4:6]);
   end
   if (inod ~= 1)
      dFdq(ic6+[4:6],ic6+[4:6]) = dFdq(ic6+[4:6],ic6+[4:6]) ...
          + TF0F*[zeros(1,3);zeros(1,3);                    ...
                  dMdthz*(dthdx*dxdthn + dthdy*dydthn)      ...
          +       dMdthzd*(dthddx*dxdthn + dthddy*dydthn)   ...
          +       dMdthzd*(dthddxd*dxddthn + dthddyd*dyddthn)];
   end
   dFdqd(ic6+[4:6],4:6) = dFdqd(ic6+[4:6],4:6)              ...
          + TF0F*[zeros(1,3);zeros(1,3);                    ...
                  dMdthzd*(dthddxd*dxddthBd + dthddyd*dyddthBd)];
   if (inod ~= 1)
      dFdqd(ic6+[4:6],ic6+[4:6]) = dFdqd(ic6+[4:6],ic6+[4:6]) ...
          + TF0F*[zeros(1,3);zeros(1,3);                      ...
                  dMdthzd*(dthddxd*dxddthnd + dthddyd*dyddthnd)];
   end

end
