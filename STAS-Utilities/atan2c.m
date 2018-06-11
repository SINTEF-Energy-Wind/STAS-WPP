function phi = atan2c (y,x);
% Works with complex step.  Protects against (0,0). The logical
% branches follow the real part of x and y.

N = size(y,1);
phi = zeros(N,1);
ry = real(y);
rx = real(x);
%del = sqrt(eps);
%ind = (abs(ry) > del) & (abs(rx) < del);  % (Inf handled OK by atan.)
%phi(ind) = sign(y(ind))*pi/2;
ind = (ry ~= 0) | (rx ~= 0);
phi(ind) = atan(y(ind)./x(ind));
inr = (rx < 0);
i2 = (inr & (ry >= 0));
i3 = (inr & (ry < 0));
phi(i2) = phi(i2) + pi;
phi(i3) = phi(i3) - pi;