function r = distributer (N,Ri,Ro,troot,ttip)

%troot = -0.35*pi;
%ttip = 0.45*pi;

t = ([0:N-1].')*(ttip - troot)/(N-1) + troot;

% The parameter t is the angle in radians.  Nodes are spaced
% according to a sine function.
r = Ri + (Ro - Ri)*(sin(t) - sin(troot))./(sin(ttip) - sin(troot));
