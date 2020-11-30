function [Vstar,phi] = clusdt (phim,Vdat,Lam,mus,sig,bet)
%
% Update the cluster wind speed estimate as well as the probability
% distribution.
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
% phim            : Vector, length Nmu (number of bins).  Probability
%                   distribution of cluster mean.
% Vdat            : Vector, length Nt.  Present wind speed at each
%                   turbine in the cluster.
% Lam             : Transition probability, phi(k+1) = Lam phi(k).
% mus             : Bins for the cluster mean wind speed.
% sig             : Width of the central part of the phi(V|mu) profile.
% bet             : Background probability of the phi(V|mu) profile.
%
% Outputs:
% --------
% Vstar           : Updated wind speed estimate at the central turbine.
% phi             : Updated probability distribution over mu's.

% First update the estimate.
phi1 = Lam*phim;

% Then incorporate the measurements.
phi = phimugV (Vdat,mus,sig,bet,phi1);

ind = (phi < 1e-12);
phi(ind) = 0.0;
phi = phi/sum(phi);  % Normalize to prevent drift from sum(p) = 1.

% Estimate of cluster wind speed.
Vstar = sum(mus.*phi);  % Expected value.

