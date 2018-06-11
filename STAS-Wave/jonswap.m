function Swf = jonswap (Hs,Tp,gamma,sigma,f)

% Pierson-Moskowitz.
Spm = ochiHubble (Hs,Tp,1,f);

Swf = (1 - 0.287*log(gamma)) ...
   .* Spm                    ...
   .* gamma.^(exp(-((f*Tp - 1)./(2*sigma)).^2));