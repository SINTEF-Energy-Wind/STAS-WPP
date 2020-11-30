function phi = phimugV (Vm,mus,sig,bet,phim)
%
% Computes phi(mu|Vm) = phi(Vm|mu) phi(mu) / Sum_k { phi(Vm|mu_k) phi(mu_k) }
%
% Version:        Changes:
% --------        -------------
% 30.06.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 30.06.2020      
%
% Inputs:
% -------
% Vm              : Set of measured values.
% mus             : Possible values of mu, must encompass range of Vm's.
% sig             : Standard deviation of the Gaussian part.
% bet             : Constant "background" probability density.
% phim            : Prior phi(mu), vector of length Nmu.
%
% Outputs:
% --------
% phi             : = phi(mu|Vm), a vector of length Nmu.

phiV = phiVgmu (Vm,mus,sig,bet);

% Normalize using mval for numerical reasons.
mval = max(max(phiV));

pphi = prod(phiV/mval) .* (phim.');
phi = (pphi.')/sum(pphi);
