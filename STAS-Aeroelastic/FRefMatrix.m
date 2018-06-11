function [F,TB_B0,dTB_B0,TB0_g] = FRefMatrix (th,th0)
% Version:        Changes:
% --------        -------------
% 08.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 08.10.2017      Double-checked formulas against report.

TB0_g = TFromTheta (th0);
[TB_B0,dTB_B0] = dTdth (th);
TB0_B = TB_B0.';
Tg_B0 = TB0_g.';
TT = TB0_B*Tg_B0;
F = [TB0_g*dTB_B0(:,1:3)*TT ...
     TB0_g*dTB_B0(:,4:6)*TT ...
     TB0_g*dTB_B0(:,7:9)*TT];

