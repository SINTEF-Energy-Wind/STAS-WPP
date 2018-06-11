function [Cld,C,Dy] = ClCdCm (Co,aoa,aq,az)
%
% Dynamic lift coefficient and standard drag and moment coefficients.
% Don't compute the coefficients from within this function, since
% this involves a lookup including a set of airfoil tables, and there
% is no need to require these tables as input to the present function.
%
%   States:           y vector:         u vector:
%   aoa               aq
%
% Version:        Changes:
% --------        -------------
% 21.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 21.12.2017      Checked using finite difference, Cld output.
%
% Inputs:
% -------
% Co              : Cl at aoa, Cd and Cm at aq, dCl/da, dCd/da, dCm/da
%                   evaluated at these points.
% aoa,aq          : Dynamic and quasi-steady angles-of-attack.
% az              : Zero-lift aoa for dynamic stall.
%
% Outputs:
% --------
% Cld             : Dynamic Cl.
% C,Dy            : d[Cl,Cd,Cm] = C*aoa + Dy*aq.

C  = zeros(3,1);
Dy = zeros(3,1);

aqa = aq - aoa;
aaz  = aoa - az;
aqaz1 = 1 + aqa/aaz;

if (absc(aaz) < eps^(0.25))
   Cld = Co(4)*aqa;
   C(1) = -Co(4);
   Dy(1) = Co(4);
else
   Cld = Co(1)*aqaz1;
   C(1) = Co(4)*aqaz1 - Co(1)/aaz - Co(1)*aqa/(aaz^2);
   Dy(1) = Co(1)/aaz;
end

%{
'---Cl---'
aqa
aaz
Co(1)
Cld
%}

Dy(2) = Co(5);
Dy(3) = Co(6);
