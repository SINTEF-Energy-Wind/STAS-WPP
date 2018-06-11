function [F,Tnn0,dTnn0] = Fmatrix (th,TBg,Tn0B)
% Caution, the input TBg = T_B0^g*T_B^B0.
% In the report this is labelled as the G matrix.
%
% Version:        Changes:
% --------        -------------
% 26.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 26.10.2017      Double-checked formulas against equations in report.
%

[Tnn0,dTnn0] = dTdth (th);
T2 = TBg*Tn0B;
F = [T2*dTnn0(:,1:3)*((T2*Tnn0).') ...
     T2*dTnn0(:,4:6)*((T2*Tnn0).') ...
     T2*dTnn0(:,7:9)*((T2*Tnn0).')];

