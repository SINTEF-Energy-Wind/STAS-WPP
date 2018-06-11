function [psie,Uzts,Dy] = rotorUZTS (xer,Ur)
%
% Axial, tangential, and spanwise components of rotorplane velocity.
% This is the local flow velocity, including the rotor's rotation.
% Can also be used for Fxyz-to-Fzts.
%
%   States:           y vector:         u vector:
%                     xer   1:3
%                     Ur    4:6
%
% Version:        Changes:
% --------        -------------
% 27.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 27.12.2017      Checked with complex step using Uzts output.
%
% Inputs:
% -------
% xer             : Position of element centroid in rotor coordinates.
% Ur              : Velocity in xyz rotor coordinates.
%
% Outputs:
% --------
% psie            : Azimuth angle of the element in rotor coordinates.
% Uzts            : Velocity in zts rotor coordinates.
% Dy              : State-space matrix.

Uzts  = zeros(3,1);
Dy    = zeros(3,6);

xmag  =  sqrt(xer(1)^2 + xer(2)^2);
xm3   =  xmag^3;
cp    =  xer(1)/xmag;
sp    =  xer(2)/xmag;
dcpdx = (xer(2)^2)/xm3;
dcpdy = -xer(1)*xer(2)/xm3;
dspdx =  dcpdy;
dspdy = (xer(1)^2)/xm3;

psie = atan2c (xer(2),xer(1));

Uzts(1) =  Ur(3);
Uzts(2) = -Ur(1)*sp + Ur(2)*cp;
Uzts(3) =  Ur(1)*cp + Ur(2)*sp;

Dy(:,1) = [0;-Ur(1)*dspdx+Ur(2)*dcpdx;Ur(1)*dcpdx+Ur(2)*dspdx];
Dy(:,2) = [0;-Ur(1)*dspdy+Ur(2)*dcpdy;Ur(1)*dcpdy+Ur(2)*dspdy];
%Dy(:,3) = 0;
Dy(:,4) = [0;-sp;cp];
Dy(:,5) = [0;cp;sp];
Dy(:,6) = [1;0;0];

%{
'---Uzts---'
psie
cp
sp
Ur
Uzts
'----------'
%}