function f = prandtl (Uz,Ut,r,Dia)
%
% Version:        Changes:
% --------        -------------
% 11.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 11.12.2017      
%

NN = size(Uz,1);

phi = atan2c(Uz,Ut);
sp = sin(phi);

f = zeros(NN,1);
ind = (real(sp) < 1e-3 & real(sp) >= 0) | (real(sp) < 0);
nind = ~ind;
snind = sum(nind);
f(ind) = 1;
f(nind) = maxc((2/pi)*acos(exp(-3*(0.5*Dia(nind) - r(nind)) ...
        ./(2*r(nind).*sp(nind)))),0.1*ones(snind,1));




