function C = dampingCNL (ces,dmu)
%
% Build the damping matrix if the damping is to be accounted for in
% the full model, as opposed to specifying an effective modal damping
% ratio on the mode-reduced model.
%
% Version:        Changes:
% --------        -------------
% 06.03.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 06.03.2018      Rotating cantilever case gives reasonable results.
%                 Same with wind turbine analysis.
%
% Inputs:
% -------
% ces             : Element damping matrix.
% dmu             : dmu/dq.
%
% Outputs:
% --------
% C               : 18-by-18 damping matrix.

C = zeros(18,18);

% Say that no damping is associated with the body reference DOFs.
% Let it all be elastic.
C(7:18,7:18) = (dmu.')*ces*dmu;

% (Same as above.)
%for kk = 1:12
%   dmkc = (dmu(:,kk).')*ces;
%   for jj = 1:12
%      C(kk+6,jj+6) = dmkc*dmu(:,jj);
%   end
%end
