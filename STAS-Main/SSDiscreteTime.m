function [A,B,dt] = SSDiscreteTime (meth,aa,bb,dti)
%
% Creates a set of discrete-time state matrices based on continuous-time
% input matrices.
%
% Version:        Changes:
% --------        -------------
% 06.10.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 06.10.2020      
%
% Inputs:
% -------
% meth            : 1 = RK2.
%                   2 = implicit trapezoidal.
% aa,bb,cc,dd     : Continuous-time state matrices.
% dti             : Prescribed timestep.  Set this to a negative number
%                   if you want the timestep generated automatically.
%
% Outputs:
% --------
% A, B, C, D      : Discrete-time state matrices.
% dt              : Timestep.

if (dti <= 0)
   % Generate the maximum frequency and timestep.
   [slap,shp,ifrq] = eigVal_silent (aa);
   fmax = abs(imag(slap(1)))/(2*pi);
   dt = 0.125/fmax;
else
   dt = dti;
end

% Generate discrete-time matrices by applying the selected integration
% method.
Nx = size(aa,1);
II = speye(Nx);
aadt = aa*dt;
bbdt = bb*dt;
if (meth == 1)
   A = II + (II + 0.5*aadt)*aadt;
   B = (II + 0.5*aadt)*bbdt;
elseif (meth == 2)
   IAp = II + 0.5*aadt;
   IAm = II - 0.5*aadt;
   A = IAm\IAp;
   B = IAm\bbdt;
end

