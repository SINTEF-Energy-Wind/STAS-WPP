function [Tn_y,Th_d,Tb_h] = basicTransforms (delta,phi)
% This builds the basic 3-by-3 transform matrices,
% associated with changes in the orientation of the
% elements in a body, in particular the nacelle turret
% and hub/pitch bearings.
%
% Version:        Changes:
% --------        -------------
% 05.12.2014      Original code.
% 30.06.2017      Adapted for complex step.
%
% Version:        Verification:
% --------        -------------
% 05.12.2014      Double-checked all the transforms.
% 30.06.2017      
%
% Inputs:
% -------
% delta       : shaft tilt angle.
% phi         : cone angle.
%
% Outputs:
% --------
% Coordinate transforms, 3-by-3, except for Th_d,
% which is 3-by-9, with a 3-by-3 block for each
% blade.

Th_d = zeros(3,9);

cd = cos(delta);
sd = sin(delta);
ch2 = cos(2*pi/3);
sh2 = sin(2*pi/3);
ch3 = cos(4*pi/3);
sh3 = sin(4*pi/3);
cphi = cos(phi);
sphi = sin(phi);

Tn_y = [ 0 sd  cd; ...
         1 0   0 ; ...
         0 cd -sd];
Th_d(:,1:3) = eye(3);
Th_d(:,4:6) = [ch2 -sh2 0; ...
               sh2  ch2 0; ...
                0    0  1];
Th_d(:,7:9) = [ch3 -sh3 0; ...
               sh3  ch3 0; ...
                0    0  1];
% Positive phi rotates the blade upwind, about the
% positive Y^b axis.
Tb_h = [cphi 0 sphi; ...
         0   1  0  ; ...
       -sphi 0 cphi];


