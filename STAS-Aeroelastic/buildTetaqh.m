function [qh,qhd,qhdd,T] = buildTetaqh (eta,etad,etadd,shape)
%
% Builds the transform from body mode to qh (constrained nodal
% displacements) degrees-of-freedom.   
%
% Version:        Changes:
% --------        -------------
% 14.01.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 14.01.2020      
%

Nr = size(shape,1);
Nn = size(shape,2);
T = [shape, sparse(Nr,Nn); sparse(Nr,Nn), shape];

qh = shape*eta;
qhd = shape*etad;
qhdd = shape*etadd;
