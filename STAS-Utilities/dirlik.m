function [rho,n] = dirlik (m0,m1,m2,m4,T,s,ds)
% This function begins with spectral moments, one set of m0, m1, m2,
% and m4.  At selected stress amplitude levels, the Dirlik approach
% is used to calculate the probability density, and, using the input
% time period, the expected number of cycles at that stress amplitude.
%
% Version:        Changes:
% --------        -------------
% 16.11.2015      Original code.
%
% Version:        Verification:
% --------        -------------
% 16.11.2015      
%
% Inputs:
% -------
% m0,m1,m2,m4     : spectral moments.
% T               : period of time over which the current spectral
%                   moments apply.
% s               : stress amplitudes.
% ds              : width of each stress amplitude.
%
% Outputs:
% --------
% rho             : Probability density of each stress amplitude.
% n               : number of cycles at each stress amplitude, over
%                   the specified time period.

sqm0 = sqrt(m0);
xm = (m1./m0).*sqrt(m2./m4);
gamma = m2./sqrt(m0.*m4);
D1 = 2*(xm - gamma.^2)./(1 + gamma.^2);
R = (gamma - xm - D1.^2)./(1 - gamma - D1 + D1.^2);
D2 = (1 - gamma - D1 + D1.^2)./(1 - R);
D3 = 1 - D1 - D2;
Q = 1.25*(gamma - D3 - D2.*R)./D1;
c1 = D1./(Q.*sqm0);
c2 = D2./(sqm0.*R.^2);
c3 = D3./sqm0;
tsqm4m2 = T*sqrt(m4./m2);

Ns = size(m0,1);
rho = sparse(size(s,1),Ns);
n = sparse(size(s,1),Ns);
for is = 1:Ns
   Z = s/sqm0(is);
   rho(:,is) = c1(is)*exp(-Z/Q(is)) ...
             + c2(is)*Z.*exp(-Z.^2/(2*R(is)^2)) ...
             + c3(is)*Z.*exp(-Z.^2/2);
   n(:,is) = rho(:,is)*ds(is)*tsqm4m2(is);
end

