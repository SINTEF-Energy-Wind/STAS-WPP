function Sv = vSpectra (vcas,dt)
%
% Perform an update of the spectra from the v cascade.
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
% vcas            : The turbulent wind speed cascade.
%
% Outputs:
% --------
% SV              : The spectra.
%

NvC = size(vcas,1);
mir = [vcas;vcas(NvC:-1:1,:)];  % Mirror to make it periodic.
Nf = 2*NvC;
df = 1/(Nf*dt);

Fv = fft(mir)/Nf;
Sv = Fv.*conj(Fv)/df;

