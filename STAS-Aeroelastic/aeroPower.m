function [Pow,y,Du] = aeroPower (qpsi,dqpsidt,Fpsi,P,idofs,blxdof,dret,Nbnod)
%
% Compute the overall aerodynamic power (force-times-velocity), as
% well as the global force and velocity associated with each blade
% node. 
%
%   States:           y vector:         u vector:
%                     vgpsi    Ndj      qpsi        Ndj
%                     Fgpsi    Ndj      dqpsi/dt    Ndj
%                                       Fpsi        Ndj
%
% Version:        Changes:
% --------        -------------
% 20.11.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 20.11.2019      Derivatives verified by complex step.  Power verified.
%
% Inputs:
% -------
% 
%
% Outputs:
% --------
% P               : Mean aero power.
% y               : v_g^psi and F_g^psi.  Power is
%                   P = 3     (v_g^psi)_0^T (F_g^psi)_0^T
%                     + (3/2) (v_g^psi)_c^T (F_g^psi)_c^T
%                     + (3/2) (v_g^psi)_s^T (F_g^psi)_s^T
% ddy             : = dy/du

Ndj = size(qpsi,1);
Nu = 3*Ndj;
Ny = 2*Ndj;

Du = sparse(Ny,Nu);

[b1,b2,b3] = MBCindices_Ndj (Ndj,idofs);
[TpsiB_Ndj,TBpsi_Ndj] = MBC (Ndj,b1,b2,b3,0);

% To get global blade node velocities in MBC coordinates, transform
% to blade coordinates, build the D matrix, then transform back.
qB = TpsiB_Ndj*qpsi;
dqBdt = TpsiB_Ndj*dqpsidt;

FB = TpsiB_Ndj*Fpsi;

vg0 = zeros(Ndj,1);
Fg0 = zeros(Ndj,1);
Dvg = sparse(Ndj,2*Ndj);   % 2*Ndj: q, dq/dt.
DFg = sparse(Ndj,2*Ndj);   % 2*Ndj: q, FB.
for inod = 2:Nbnod
   
   for ib = 1:3

      rdof = idofs(5+ib);
      dof  = idofs(5+ib) + 6*(inod-1);

      % Get the global velocity in global coordinates.
      qr = qB(rdof+[1:6]);
      qn = qB(dof+[1:6]);
      Pr = P(rdof+[1:6]);
      Pn = P(dof+[1:6]);
      dqrdt = dqBdt(rdof+[1:6]);
      dqndt = dqBdt(dof+[1:6]);
      [vg,ddy] = globalVelocity (qn,qr,Pn,Pr,dqndt,dqrdt);
      vg0(dof+[1:6]) = vg;
      Dvg(dof+[1:6],[rdof+[1:6] dof+[1:6] Ndj+rdof+[1:6] Ndj+dof+[1:6]]) = ddy;

      % Express the blade forces in global coordinates.
      [TBB0,dTBB0] = dTdth (qr(4:6));
      TB0g = TFromTheta (Pr(4:6));
      T = TB0g*TBB0;
      dT = TB0g*dTBB0;
      TT = [T, zeros(3,3);zeros(3,3), T];
      Fg0(dof+[1:6]) = TT*FB(dof+[1:6]);
      
      DFg(dof+[1:6],Ndj+dof+[1:6]) = TT;
      for jj = 1:3
         jc3 = 3*(jj-1);
         dTT = [dT(:,jc3+[1:3]), zeros(3,3);zeros(3,3), dT(:,jc3+[1:3])];
         DFg(dof+[1:6],rdof+3+jj) = dTT*FB(dof+[1:6]);
      end

%'---------------------'
%inod
%ib
%[qr qn Pr Pn dqrdt dqndt vg]
%ddy
%[FB(dof+[1:6]) Fg0(dof+[1:6])]
%T

   end

end

% Transform back to MBC.
TpsiB_2Ndj = [TpsiB_Ndj, sparse(Ndj,Ndj); sparse(Ndj,Ndj), TpsiB_Ndj];
y = [TBpsi_Ndj*vg0;TBpsi_Ndj*Fg0];
Du(1:Ndj,1:2*Ndj) = TBpsi_Ndj*Dvg*TpsiB_2Ndj;
Du(Ndj+[1:Ndj],[[1:Ndj] 2*Ndj+[1:Ndj]]) = TBpsi_Ndj*DFg*TpsiB_2Ndj;

Pow = sum(Fg0.*vg0);