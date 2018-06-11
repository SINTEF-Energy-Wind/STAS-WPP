function [Fr,Dy] = forcesRotorCS (Fa,Tar,dTar)
%
% The airfoil forces in rotorplane coordinates.
%
%   States:           y vector:         u vector:
%                     qy   6
%                     qp   6
%                     qn1  6
%                     qn2  6
%                     Fa   6
%
% Version:        Changes:
% --------        -------------
% 02.01.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 02.01.2017      Verified using complex step and Fr output.
%
% Inputs:
% -------
% Tar,dTar        : Transform from airfoil to rotor: qy,qp,qn1,qn2.
% Fa              : Local forces in airfoil coordinates.
%
% Outputs:
% --------
% Fr              : Local forces in rotor coordinates.
% Dy              : State space.

Dy = zeros (6,30);

mat = [Tar zeros(3,3);zeros(3,3) Tar];
Fr = mat*Fa;

Dy(:,25:30) = mat;

for jj = 1:24
   jc3 = 3*(jj-1);
   Dy(:,jj) = [dTar(:,jc3+[1:3]) zeros(3,3);zeros(3,3) dTar(:,jc3+[1:3])]*Fa;
end
