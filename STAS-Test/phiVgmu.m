function phi = phiVgmu (Vm,mus,sig,bet)
%
% Compute p(V_m|mu_V), the probability of a set of measured values Vm
% given a set of possible values mus.
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
%
% Outputs:
% --------
% phi             : p(V_m|mu_V), NV-by-Nmu

NV = size(Vm,1);
Nmu = size(mus,1);

del = Vm - mus.';
%delr = [[mus(1)-mus(Nmu:-1:2)];[mus-mus(1)]];  % Range of del.
%Ndel = size(delr,1);
%w = delr(Ndel) + 0.5*(delr(Ndel)-delr(Ndel-1)) ...
%  - delr(1)    - 0.5*(delr(2)-delr(1));
w = mus(Nmu) + 0.5*(mus(Nmu)-mus(Nmu-1)) ...
  - mus(1)   - 0.5*(mus(2)-mus(1));

%pbet = w*bet;  % Total probability of the background part.
%phi = bet + (1 - pbet)*normpdf(del,0.0,sig);

phi = (bet/w) + (1 - bet)*normpdf(del,0.0,sig);


