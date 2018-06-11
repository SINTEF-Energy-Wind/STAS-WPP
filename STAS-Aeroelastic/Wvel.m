function [W,C,Dy] = Wvel (Viz,Vzts,f)
%
% Linearization of the local rotorplane velocity magnitude W.
%
%   States:           y vector:         u vector:
%   Viz               Vzts   1:3
%                     f       4
%
% Version:        Changes:
% --------        -------------
% 02.01.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 02.01.2017      Verified by complex step using W output.
%
% Inputs:
% -------
% Viz             : Z-component induced velocity (state).
% Vzts            : Local rotorplane flow velocity, zts components.
% f               : Prandtl factor.
%
% Outputs:
% --------
% W               : Local rotorplane velocity magnitude.
% C,Dy            : State matrices.

Dy = zeros(1,4);

vfv = Vzts(1) + f*Viz;
W = sqrt(vfv^2 + Vzts(2)^2 + Vzts(3)^2);

C     = vfv*f/W;
Dy(1) = vfv/W;
Dy(2) = Vzts(2)/W;
Dy(3) = Vzts(3)/W;
Dy(4) = vfv*Viz/W;

