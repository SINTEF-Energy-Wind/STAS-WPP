function [zeta0,eps,k,df] = waveFourier (depth,Hs,Tp,f)
%
% Computes the amplitude and randomized phase for each of a number
% of specified frequencies.  Also calls a subroutine to calculate
% k from the dispersion relationship.
%
% Version:        Changes:
% --------        -------------
% 12.03.2015      Adapted from waveFourier.f90.
%
% Version:        Verification:
% --------        -------------
% 12.03.2015      Memo AN 15.12.19
%
% Inputs:
% -------
% depth           : Water depth.
% Hs              : Significant wave height.
% Tp              : Wave period at peak energy.
% f               : List of frequencies.
%
% Outputs:
% --------
% zeta0           : Array of wave amplitudes.
% eps             : Array of random phases.
% k               : Wave number.
% df              : Frequency bin widths.

Nf = size(f,1);

df = getdf(f);

gamma = DNVgamma(Hs,Tp);

sigma = DNVsigma(f,Tp);

% Jonswap:
Szeta = jonswap (Hs,Tp,gamma,sigma,f);

% Pierson-Moskowitz:
%Szeta = ochiHubble (Hs,Tp,1,f);

% Example of a dual-peak spectrum.
%Szeta = ochiHubble (0.97,6.4,0.8,f) + ochiHubble (0.72,16,0.3,f);

% Clear numerical noise in the sparse spectrum.
Szeta(abs(Szeta)<1e-12) = 0;

eps = 2.d0*pi*rand(Nf,1);

zeta0 = sqrt(2*Szeta.*df);

k = getk(f,depth);



% Possible output of interest.
%fid = fopen('spec.txt','w');
%for ifreq = 1:Nf
%   fprintf(fid,'%+5.4e %+5.4e\n',f(ifreq),Szeta(ifreq));
%end
%fclose(fid);
