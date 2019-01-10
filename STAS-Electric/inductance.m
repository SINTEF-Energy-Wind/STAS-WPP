function [L,dL] = inductance (I0,L0,II)
%
% Saturation effects in magnetic components are represented by a 
% nonlinear effective inductance, that at any given point is the tangent
% of the terminal voltage versus rate-of-change-of-current curve.  This
% curve is represented by a spline.
%
% Version:        Changes:
% --------        -------------
% 27.08.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 27.08.2018      Checked sample test cases.
%
% Inputs:
% -------
% I0              : I's at which L's are given, four points.
% L0              : values of inductance at the I's.
% II              : currents at which to compute L.
%
% Outputs:
% --------
% L               : effective value of inductance.
% dL              : = dL/dI.

L = zeros(size(II,1),1);
dL = zeros(size(II,1),1);

[A,B] = splineCoefficients (I0);
[S,dSdx] = splineSMatrix (I0,II);
Am1BL0 = A\(B*L0);
L = S*Am1BL0;
dL = dSdx*Am1BL0;

