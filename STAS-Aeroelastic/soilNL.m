function F = soilNL (kx,ky,kz,kthz,cx,cy,cz,cthz,PF,qF,dqFdt)
%
% Add soil stiffness and damping to the foundation.  Consider the
% stiffness.  Each nodal spring applies a force F^f0 = -ks*d^f0,
% with d^f0 the lateral or vertical displacement in undeformed
% foundation coordinates.  Together, the stiffness and damping
% nodal forces can be written
% F^f0 = T_F0^F*K*(T_g^F0*O_F/F0^g + T_F^F0(P_n0/F^F + q_n/n0^F) - P_n0/F^F)
%      + T_F0^F*C*(T_g^F0*dO_F/F0^g/dt + T_F^F0*dq_n/n0^F/dt
%      +           dT_F^F0/dPhi*(Pn0/F^F + q_n/n0^F)*dPhi/dt).
% The torsional displacement is defined as 
% thz = atan(y/x), where x = (TF^F0*Tn0^F*T_n^n0)(1,2) and z = ()(2,2),
% with torsional moment Mz^F = -kz*thz.
% The torsional velocity is dthz/dt, with torsional moment
% Mz^F = 
% Applying the nodal forces requires the operation Qh*F.  This
% is accomplished outside this function, which returns the nodal
% force.
%
% Version:        Changes:
% --------        -------------
% 13.05.2017      Original code, reformulation of soil.m.
%
% Version:        Verification:
% --------        -------------
% 13.05.2017      Verified results for some simple cases.
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
%

Nnod = size(kx,1);
Nd = 6*Nnod;

F = zeros(Nd,1);

TF0g = TFromTheta (PF(4:6));
TgF0 = TF0g.';
[TFF0,dTFF0] = dTdth (qF(4:6));
TF0F = TFF0.';

for inod = 1:Nnod

   ic6 = 6*(inod-1);

   Kd = diag([kx(inod) ky(inod) kz(inod)].');
   Cd = diag([cx(inod) cy(inod) cz(inod)].');

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
   F(ic6+[1:3]) = -TF0F*Kd*(TgF0*qB(1:3)                        ...
                +           TFF0*(Pn(1:3) + qn(1:3)) - Pn(1:3)) ...
                -  TF0F*Cd*(TgF0*dqB(1:3) + TFF0*dqn(1:3));

   for jj = 1:3
      jc3 = 3*(jj-1);
      F(ic6+[1:3]) = F(ic6+[1:3])               ...
                   - TF0F*Cd*dTFF0(:,jc3+[1:3]) ...
                   * (Pn(1:3) + qn(1:3))*dqB(3+jj);
   end

   % Moments due to torsion.
   Tn0F = TFromTheta (Pn(4:6));
   [Tnn0,dTnn0] = dTdth (qn(4:6));
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

end
