function [Vzts,Dy] = rotorVZTS (Uzts,Vir,r,Omega)
%
% Axial, tangential, and spanwise components of rotorplane velocity.
% This is the local flow velocity, not including the rotor's rotation.
%
%   States:           y vector:         u vector:
%                     Uzts  1:3
%                     Vir   4:5
%                     r      6
%                     Omega  7
%
% Version:        Changes:
% --------        -------------
% 27.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 27.12.2017      Checked some cases, Vzts matches the equivalent
%                 value in BEMNL.m.
%
% Inputs:
% -------
% 
%
% Outputs:
% --------
% 

Vzts = zeros(3,1);
Dy   = zeros(3,7);

Vzts(1) = Uzts(1) - Vir(1);
Vzts(2) = Uzts(2) - Vir(2) + r*Omega;
Vzts(3) = Uzts(3);

Dy(:,1:3)   =  eye(3);
Dy(1:2,4:5) = -eye(2);
Dy(2,6)     =  Omega;
Dy(2,7)     =  r;

%{
'---Vzts---'
Uzts
Vir
r*Omega
Vzts
'--end Vzts--'
%}
