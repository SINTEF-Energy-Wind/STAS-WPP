function [rho,n] = rayleigh (m0,m2,T,s,ds)
% This function begins with spectral moments, one set of m0, m1, m2,
% and m4.  At selected stress amplitude levels, the Rayleigh formula
% is used to calculate the probability density, and, using the input
% time period, the expected number of cycles at that stress amplitude.
%
% Version:        Changes:
% --------        -------------
% 22.09.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 22.09.2020      
%
% Inputs:
% -------
% m0,m2           : spectral moments.
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

sig2 = m0;
fc = sqrt(m2/m0);
rho = (s./sig2).*exp(-0.5*(s.^2)/sig2);
n = rho.*ds.*fc*T;



