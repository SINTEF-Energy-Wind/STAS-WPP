function [Fld,Fa,Dy] = airfoilForces (rho,ch,L,x,y,aq,U,Cl,Cd,Cm)
%
% Compute the airfoil lift, drag, and moment.
%
%   States:           y vector:         u vector:
%                     aq
%                     U
%                     Cl,Cd,Cm
%
% Version:        Changes:
% --------        -------------
% 21.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 21.12.2017      Checked by complex step using Fa.
%
% Inputs:
% -------
% rho             : Air density.
% ch              : Chord.
% L               : Spanwise length.
% x,y             : Location of the origin of the section coordinate system
%                   with respect to the point of definition of the moment
%                   coefficient, typically the quarter-chord.
%
% Outputs:
% --------
% 

pclu = 0.5*rho*ch*L*(U^2);
Fld = pclu*[Cl;Cd;ch*Cm];

sa = sin(aq);
ca = cos(aq);

Tldm = [1 0 0;0 1 0;0 0 ch];
T = [-sa ca 0;ca sa 0;0 0 0;0 0 0;0 0 0;-ca*x-sa*y -sa*x+ca*y -ch];
dTda = [-ca -sa 0; -sa ca 0; 0 0 0; 0 0 0; 0 0 0; sa*x-ca*y -ca*x-sa*y 0];

C = [Cl;Cd;Cm];
Fa = pclu*T*C;

Dy = zeros(9,5);
Dy(1:3,2)   = 0.5*rho*ch*L*2*U*Tldm*C;
Dy(1:3,3:5) = pclu*Tldm;
Dy(4:9,1)   = pclu*dTda*C;
Dy(4:9,2)   = 0.5*rho*ch*L*2*U*T*C;
Dy(4:9,3:5) = pclu*T;
