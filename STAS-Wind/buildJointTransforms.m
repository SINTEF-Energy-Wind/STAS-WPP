function [Ty_g,Td_n,Tp_b] = buildJointTransforms (chi,psi,beta)
% This builds the basic 3-by-3 transform matrices,
% associated with the joint degrees-of-freedom.  This
% does not include elastic deformation.
%
% Version:        Changes:
% --------        -------------
% 05.12.2014      Original code.
% 02.08.2017      Adapted for complex step derivatives.
%
% Version:        Verification:
% --------        -------------
% 05.12.2014      Double-checked all the transforms.
% 02.08.2017      
%
% Inputs:
% -------
% chi         : yaw angle.
% psi         : azimuth angle.
% beta        : vector length 3, pitch angles.
%
% Outputs:
% --------
% Coordinate transforms, 3-by-3, except for Tp_b, which
% is 9-by-3, with a 3-by-3 block for each blade.

Tp_b = zeros(9,3);

cchi = cos(chi);
schi = sin(chi);
cp   = cos(psi);
sp   = sin(psi);
cb1  = cos(beta(1));
sb1  = sin(beta(1));
cb2  = cos(beta(2));
sb2  = sin(beta(2));
cb3  = cos(beta(3));
sb3  = sin(beta(3));

Ty_g = [cchi -schi 0; ...
        schi  cchi 0; ...
         0     0   1];
Td_n = [cp -sp 0; ...
        sp  cp 0; ...
        0   0  1];
% Positive beta is rotation about the NEGATIVE X^b axis.
Tp_b(1:3,:) = [1   0   0;  ...
               0  cb1 sb1; ...
               0 -sb1 cb1];
Tp_b(4:6,:) = [1   0   0;  ...
               0  cb2 sb2; ...
               0 -sb2 cb2];
Tp_b(7:9,:) = [1   0   0;  ...
               0  cb3 sb3; ...
               0 -sb3 cb3];

