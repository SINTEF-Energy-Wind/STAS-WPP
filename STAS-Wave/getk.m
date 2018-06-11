function k = getk (f,d)
%
% This returns the value of k that satisfies k tanh kd = w^2/g.
% The function is the leading behavior sqrt(x) plus a polynomial
% fit through x=5, above which y=x is a very close approximation.
%
% Inputs:
% -------
% f               : Frequencies in Hz.
% d               : Water depth.
%
% Outputs:
% --------
% k               : Wave number.

w = 2*pi*f;
x = d*(w.^2)/9.807;

k = (sqrt(x) ...
  + ((((0.000380314259937.*x ...
  -     0.00418571290767).*x ...
  +     0.00275998776826).*x ...
  +     0.142445384227  ).*x ...
  +     0.0573461092478 ).*x)/d;

k(x >= 5) = x(x >= 5)/d;