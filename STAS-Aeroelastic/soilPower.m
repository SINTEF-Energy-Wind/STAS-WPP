function [y,Du] = soilPower (s,qf,dqfdt,P)
%
% Compute the power transmitted through the soil springs. 
%
%   States:           y vector:         u vector:
%                     vg       6*Nnod   qf        6*Nnod
%                     Fspr     6*Nnod   dqfdt     6*Nnod
%
% Version:        Changes:
% --------        -------------
% 20.11.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 20.11.2019      Derivatives verified by complex step.
%
% Inputs:
% -------
% s               : Data structure from STASTurbine input.
% qf,dqfdt        : Nodal displacements and velocities of the
%                   foundation.
% P               : Vector of nodal offsets for the foundation.
%
% Outputs:
% --------
% y               : [vg;Fsoil]
% Du              : dy/d[q;dqdt].

Ndf = size(qf,1);
Nu = 2*Ndf;
Ny = 2*Ndf;
Nnod = s.foundation.Nnod;

% ksoil: [kxg kyg kzg kthzg cxg cyg czg cthzg].' for each node.
ksoil = s.foundation.k.';

Du = sparse(Ny,Nu);
y = zeros(Ny,1);

[F,dFdq,dFdqd] = soilLin (ksoil(:,1),ksoil(:,2),ksoil(:,3),ksoil(:,4), ...
                          ksoil(:,5),ksoil(:,6),ksoil(:,7),ksoil(:,8), ...
                          P,qf,dqfdt);

y(Ndf+[1:Ndf]) = F;
Du(Ndf+[1:Ndf],[1:Ndf]) = dFdq;
Du(Ndf+[1:Ndf],Ndf+[1:Ndf]) = dFdqd;

rdof = 0;
TB0g = TFromTheta (P(rdof+[4:6]));

% Reference node.
qn = zeros(6,1);
qB = qf(rdof+[1:6]);
Pn = zeros(6,1);
PB = P(rdof+[1:6]);
dqndt = zeros(6,1);
dqBdt = dqfdt(rdof+[1:6]);
[xng,ddyp] = globalPosLin (qB,PB,qn,Pn);
[vg,ddyv] = globalVelocity (qn,qB,Pn,PB,dqndt,dqBdt);

y(rdof+[1:6]) = vg;
Du(rdof+[1:6],[rdof+[1:6] Ndf+rdof+[1:6]]) = ddyv(:,[[1:6] [13:18]]);

for inod = 2:s.foundation.Nmud+1

   dof = 6*(inod-1);

   qn = qf(dof+[1:6]);
   qB = qf(rdof+[1:6]);
   Pn = P(dof+[1:6]);
   PB = P(rdof+[1:6]);
   dqndt = dqfdt(dof+[1:6]);
   dqBdt = dqfdt(rdof+[1:6]);

   [xng,ddyp] = globalPosLin (qB,PB,qn,Pn);
   [vg,ddyv] = globalVelocity (qn,qB,Pn,PB,dqndt,dqBdt);

   y(dof+[1:6]) = vg;
   indv = [rdof+[1:6] dof+[1:6] Ndf+rdof+[1:6] Ndf+dof+[1:6]];
   Du(dof+[1:6],indv) = ddyv;

end
