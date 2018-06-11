function [Ua,Du,Dy] = localFlow (Tag,dTag,Tar,dTar,Vg,Vir,wa)
%
% The flow local to the airfoil, in airfoil coordinates.
%
%   States:           y vector:         u vector:
%                     qy     1:6        Vg
%                     qp     7:12
%                     qn1   13:18
%                     qn2   19:24
%                     wa    25:27
%                     Vir   28:30
%
% Version:        Changes:
% --------        -------------
% 21.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 21.12.2017      Checked using complex step via Ua output.
%
% Inputs:
% -------
% Tag,dTag        : Transform from airfoil to global: qp,qn1,qn2.
% Tar,dTar        : Transform from airfoil to rotor: qy,qp,qn1,qn2.
% Vg,Vir,wa       : Velocities.
%
% Outputs:
% --------
% Ua              : Local flow in airfoil coordinates.
% Du,Dy           : State space, dUa = Dy*[dqy,dqp,dqn1,dqn2,dwa,dVir] + Du*dVg.

Dy = zeros (3,30);

Tga = Tag.';
Tra = Tar.';

Ua =  Tga*Vg + Tra*Vir - wa;

%{
'---Ua---'
Tga
Tra
Tgr = Tar*Tga
Vg
Vir
wa
wg = Tag*wa
Ua
Ur = Tar*Ua
'--------'
%}

Du = Tga;

Dy(:,25:27) = -eye(3);
Dy(:,28:30) = Tra;

for jj = 1:6

   j3a = 3*(jj-1);

   % qy.  Only Tar.
   Dy(:,jj) = (dTar(:,j3a+[1:3]).')*Vir;

end

for jj = 1:18

   j3a = 3*(jj+6-1);
   j3b = 3*(jj-1);

   Dy(:,jj+6) = (dTag(:,j3b+[1:3]).')*Vg + (dTar(:,j3a+[1:3]).')*Vir;

end