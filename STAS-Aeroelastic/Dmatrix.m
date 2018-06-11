function D = Dmatrix (F)
% Version:        Changes:
% --------        -------------
% 08.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 08.10.2017      Double-checked against equations in report.

D = [-F(2,3) -F(2,6) -F(2,9);
      F(1,3)  F(1,6)  F(1,9);
     -F(1,2) -F(1,5) -F(1,8)];
