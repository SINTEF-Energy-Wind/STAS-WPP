function [psie,vxyz,Dy] = ztsToxyz (xer,vzts)
%
% Transform a vector from zts to xyz coordinates.
%
%   States:           y vector:         u vector:
%                     xer   1:3
%                     vzts  4:6
%
% Version:        Changes:
% --------        -------------
% 29.01.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 29.01.2018      Verified by complex step on vxyz output.
%
% Inputs:
% -------
% xer             : Position of element centroid in rotor coordinates.
% vzts            : Vector in zts rotor coordinates.
%
% Outputs:
% --------
% psie            : Azimuth angle of the element in rotor coordinates.
% vxyz            : Vector in zts rotor coordinates.
% Dy              : State-space matrix.

vxyz = zeros(3,1);
Dy   = zeros(3,6);

xmag  =  sqrt(xer(1)^2 + xer(2)^2);
xm3   =  xmag^3;
cp    =  xer(1)/xmag;
sp    =  xer(2)/xmag;
dcpdx = (xer(2)^2)/xm3;
dcpdy = -xer(1)*xer(2)/xm3;
dspdx =  dcpdy;
dspdy = (xer(1)^2)/xm3;

psie = atan2c (xer(2),xer(1));

vxyz(1) = -vzts(2)*sp + vzts(3)*cp;
vxyz(2) =  vzts(2)*cp + vzts(3)*sp;
vxyz(3) =  vzts(1);

Dy(:,1) = [-vzts(2)*dspdx+vzts(3)*dcpdx;vzts(2)*dcpdx+vzts(3)*dspdx;0];
Dy(:,2) = [-vzts(2)*dspdy+vzts(3)*dcpdy;vzts(2)*dcpdy+vzts(3)*dspdy;0];
%Dy(:,3) = 0;
Dy(:,4) = [0;0;1];
Dy(:,5) = [-sp;cp;0];
Dy(:,6) = [cp;sp;0];
