function [T,dT] = dTdth (th)
%
% Some effort was required for this to work with complex step.  The
% problem is that the standard formula underflows if the real part
% is zero and the complex part is very small.  A power series 
% solution is used for small inputs.  This has been tested using
% a complex step of sqrt(eps).
%
% Version:        Changes:
% --------        -------------
% 08.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 08.10.2017      Checked against an independent calculation using
%                 alternate formulae.
%
% Inputs:
% -------
% th              : theta parameters.
%
% Outputs:
% --------
% T               : 3-by-3.
% dT              : dT/dth, 3-by-9.

normth = sqrt(th.'*th);

dT = zeros(3,9);

% Compute the T matrix and spin(th).
TTH = vecToSpin(th);
T = expm (TTH);

if (abs(normth) < eps^0.25)  % Complex magnitude is intentionally included
                             % in this norm.

   normth2 = normth^2;
   normth4 = normth2^2;

   c1 = -1/3 + normth2/30;
   c2 = 1 - normth2/6 + normth4/120;
   c3 = -1/12 + normth2/180;
   c4 = 0.5 - normth2/24 + normth4/720;

   for kk = 1:3
      k3 = 3*(kk-1);
      ee = zeros(3,1);
      ee(kk) = 1;
      sigk = vecToSpin(ee);
      dT(:,k3+[1:3]) = c1*th(kk)*TTH + c2*sigk ...
                     + c3*th(kk)*TTH*TTH + c4*(sigk*TTH + TTH*sigk);   
   end

else

   Tp = T.';

   % Compute the partial derivatives of the T matrix.
   % Gallego G, Yezzi A (2015).  A compact formula for the derivative
   % of a 3-D rotation in exponential coordinates.  Journal of
   % Mathematical Imaging and Vision 51: 378-384.
   for ie = 1:3
      ee = zeros(3,1);
      ee(ie) = 1;
      term1 = th(ie)*TTH;
      term2 = vecToSpin(TTH*(eye(3) - T)*ee);
      dT(:,3*(ie-1)+[1:3]) = ((term1 + term2)/(normth^2))*T;
   end   

end

