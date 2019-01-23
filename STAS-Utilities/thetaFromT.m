function th = thetaFromT (T)
%
% Felippa and Haugen (2005).  A unified formulation of small-strain
% corotational finite elements: I. Theory.  Computational Methods
% in Applied Mechanics and Engineering 194: 2285-2335. 
%
% The solution is obtained in two stages.  First, quaternions are
% obtained from the T matrix.  Then, exponential map parameters
% from the quaternions.
%
% Version:        Changes:
% --------        -------------
% 04.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 04.10.2017      Verified that this inverts TFromTheta, which is
%                 based on expm(spin(th)).
%

% Compute the quaternions.
tr = trace(T);
[md,imd] = max(diag(real(T)));

p = zeros(4,1);

if (real(tr) >= real(md))
   p(1) = 0.5*sqrt(1 + tr);
   p(2) = 0.25*(T(3,2) - T(2,3))/p(1);
   p(3) = 0.25*(T(1,3) - T(3,1))/p(1);
   p(4) = 0.25*(T(2,1) - T(1,2))/p(1);
else
   jj = imd + [1:2];
   jj(jj > 3) = jj(jj > 3) - 3;
   p(imd+1) = sqrt(0.5*T(imd,imd) + 0.25*(1 - tr));
   p(1) = 0.25*(T(jj(2),jj(1)) - T(jj(1),jj(2)))/p(imd+1);
   p(jj(1)+1) = 0.25*(T(jj(1),imd) + T(imd,jj(1)))/p(imd+1);
   p(jj(2)+1) = 0.25*(T(jj(2),imd) + T(imd,jj(2)))/p(imd+1);
end

% Compute the exponential map parameters from the quaternions.
magth = 2*acos(p(1));
th = p(2:4)/(0.5*sinc(magth/(2*pi)));


