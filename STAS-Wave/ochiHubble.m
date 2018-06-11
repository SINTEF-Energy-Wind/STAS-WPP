function Swf = ochiHubble (Hs,Tp,lambda,f)
%
% Calculates one of the Ochi-Hubble spectra.
%
% Version:        Changes:
% --------        -------------
% 12.03.2015      Adapted from ochiHubble.f90.
%
% Version:        Verification:
% --------        -------------
% 12.03.2015      
%
% Inputs:
% -------
% Hs              : Significant wave height.
% Tp              : Mean period.
% lambda          : Shape parameter.
% f               : Frequency (Hz).
%
% Outputs:
% --------
% Swf             : Wave spectrum.

lp25 = lambda + 0.25;
ftp = f*Tp;
Swf = ((Hs^2)*Tp*(lp25^lambda) ...
   ./ (4*gamma(lambda)))       ...
   .* (ftp.^(-4*lp25))         ...
   .* exp(-lp25./(ftp.^4));