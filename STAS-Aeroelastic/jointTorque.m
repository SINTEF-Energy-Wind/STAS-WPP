function [Mm,Dy] = jointTorque (qm,Pm,Tjm,Mj)
%
% Apply an actuator torque at a joint, and the reaction torque at
% the adjacent master node.
%
%   States:           y vector:         u vector:
%                     qm
%                     Mj
%
% Version:        Changes:
% --------        -------------
% 27.02.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 27.02.2018      
%
% Inputs:
% ------- 
% qm,Pm           : Master node deformations and undeformed
%                   positions.
% Tjm             : Transform from joint to master node coordinates,
%                   including the joint rotation.
% Mj              : Joint moment in body reference coordinates.
%
% Outputs:
% --------
% Mm              : Joint moment reaction at master node.
% Dy              : Linearization.

Dy = zeros(3,9);

[Tmm0,dTmm0] = dTdth (qm(4:6));
Tm0B = TFromTheta (Pm(4:6));

T = Tm0B*Tmm0*Tjm;
Mm = T*Mj;

Dy(:,7:9) = T;

for jj = 1:3
   jc3 = 3*(jj-1);
   Dy(:,jj+3) = Tm0B*dTmm0(:,jc3+[1:3])*Tjm*Mj;
end

% Assumed that the joint is aligned such that the joint moment
% does not vary with the joint angle, hence we don't need to
% consider derivatives of Tjm.
