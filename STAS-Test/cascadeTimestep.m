function [D,] = cascadeTimestep (xcin,V,Ac,Bc,df)
%
% Perform an update of the V cascade and associated outputs.
%
% Version:        Changes:
% --------        -------------
% 14.10.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 14.10.2020      
%
% Inputs:
% -------
% 
%
% Outputs:
% --------
% 
%

Nxc = size(xcin,1);
xc = Ac*xcin + Bc*V
xct = [xc; xc(Nxc:-1:1)];  % Mirror to make it periodic.
Nf = 2*Nxc;

FV = fft(xc)/Nf;
SV = FV.*conj(FV)/df;

