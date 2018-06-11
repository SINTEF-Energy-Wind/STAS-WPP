function [aq,Dy] = alphaq (Ua)
%
% The quasi-steady angle-of-attack.
%
%   States:           y vector:         u vector:
%                     Ua
%
% Version:        Changes:
% --------        -------------
% 21.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 21.12.2017      Verified by complex step on aq output.  Note, the
%                 precision is lower than machine precision, roughly
%                 10^-9 instead of 10^-16.  I suspect the atan2c
%                 function.  [Verified, it is the call to atan () in
%                 atan2c that results in loss of precision.]
%
% Inputs:
% -------
% Ua              : Local flow in airfoil coordinates.
%
% Outputs:
% --------
% aq              : QS aoa.
% Dy              : daq = Dy*[dUa,x;dUa,y]

aq = atan2c (Ua(2),Ua(1));

val = Ua(1)^2 + Ua(2)^2;
Dy = [-Ua(2)/val Ua(1)/val 0];

%{
'---aq---'
Ua
aq
'--end aq--'
%}