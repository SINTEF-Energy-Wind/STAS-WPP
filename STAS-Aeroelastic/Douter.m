function [Do,Dy] = Douter (xnr)
%
% An effective Do can be found from the projected rotor coordinates of
% the outermost node on the present blade.
%
%   States:           y vector:         u vector:
%                     xnr
%
% Version:        Changes:
% --------        -------------
% 30.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 30.12.2017      Gives reasonable values.
%

Do = 2*sqrt(xnr(1)^2 + xnr(2)^2);

Dy = 4*[xnr(1) xnr(2) 0]/Do;

