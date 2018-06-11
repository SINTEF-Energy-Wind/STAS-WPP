function TB0_g = undeformedBodyToGlobal (yaw,tilt,azimuth,cone,pitch)
%
% Computes the undeformed orientation of each body's reference
% coordinate system.
%
% Version:        Changes:
% --------        -------------
% 03.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 03.10.2017      Checked some simple cases.
%
% Inputs:
% -------
% yaw ... pitch   : Input angles.  Note that positive blade pitch 
%                   rotates the leading edge into the wind -- a 
%                   rotation about the -X^p axis.
%
% Outputs:
% --------
% TB0_g
%

TB0_g = sparse(3,21);
Tp_b  =  zeros(3,9);

[Tn_y,Th_d,Tb_h] = basicTransforms (tilt,cone);

sx = sin(yaw);
cx = cos(yaw);
Ty_g = [cx -sx   0; ...
        sx  cx   0; ...
         0   0   1];

sp = sin(azimuth);
cp = cos(azimuth);
Td_n =      [cp -sp   0; ...
             sp  cp   0; ...
              0   0   1];

sb = sin(pitch(1));
cb = cos(pitch(1));
Tp_b(:,1:3) =     [1   0   0; ...
                   0  cb  sb; ...
                   0 -sb  cb];

sb = sin(pitch(2));
cb = cos(pitch(2));
Tp_b(:,4:6) =     [1   0   0; ...
                   0  cb  sb; ...
                   0 -sb  cb];

sb = sin(pitch(3));
cb = cos(pitch(3));
Tp_b(:,7:9) =     [1   0   0; ...
                   0  cb  sb; ...
                   0 -sb  cb];

TB0_g(:,1:3)   = eye(3);         % Foundation.
TB0_g(:,4:6)   = eye(3);         % Tower.
TB0_g(:,7:9)   = Ty_g;           % Nacelle.
TB0_g(:,10:12) = Ty_g*Tn_y*Td_n; % Driveshaft.
TB0_g(:,13:15) = TB0_g(:,10:12)*Th_d(:,1:3)*Tb_h*Tp_b(:,1:3); % Blade 1.
TB0_g(:,16:18) = TB0_g(:,10:12)*Th_d(:,4:6)*Tb_h*Tp_b(:,4:6); % Blade 2.
TB0_g(:,19:21) = TB0_g(:,10:12)*Th_d(:,7:9)*Tb_h*Tp_b(:,7:9); % Blade 3.
