function [Ur,Dy] = localFlowRotorCS (Tar,dTar,Ua)
%
% The flow local to the airfoil, in rotorplane coordinates, for use
% in computing the Prandtl factor.
%
%   States:           y vector:         u vector:
%                     qy
%                     qp
%                     qn1
%                     qn2
%                     Ua
%
% Version:        Changes:
% --------        -------------
% 25.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 25.12.2017      Verified using complex step, Ur output.
%
% Inputs:
% -------
% Tar,dTar        : Transform from airfoil to rotor: qy,qp,qn1,qn2.
% Ua              : Local flow velocity in airfoil coordinates.
%
% Outputs:
% --------
% Ur              : Local flow in rotor coordinates.
% Dy              : State space, dUr = Dy*[dqy,dqp,dqn1,dqn2,dUa].
%                   qy,qp,qn1,qn2 are length 6, Ua length 3.

Dy = zeros(3,27);

Ur = Tar*Ua;

for jj = 1:24
   jc3 = 3*(jj-1);
   Dy(:,jj) = dTar(:,jc3+[1:3])*Ua;
end

Dy(:,25:27) = Tar;

