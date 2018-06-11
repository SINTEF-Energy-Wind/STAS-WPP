function [dasdt,A,By] = aoadyn (ch,as,aq,ast,U)
%
% The dynamic angle-of-attack.
%
%   States:           y vector:         u vector:
%   aoa               aq
%   a1                Umag
%   a2
%
% Version:        Changes:
% --------        -------------
% 21.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 21.12.2017      Checked by complex step, using dasdt.  At present the
%                 scheduling of tau at deep-stall AOA's is deactivated.
%
% Inputs:
% -------
% ch              : Chord length.
% as              : [aoa,a1,a2].
% aq              : QS aoa.
% ast             : Deep-stall angle of attack.
% U               : Magnitude of local flow velocity.
%
% Outputs:
% --------
% dasdt           : d/dt[aoa,a1,a2]
% A               : 3-by-3 state matrix.
% By              : 3-by-2 matrix.

%{
'---------aoadyn----------'
ch
as
aq
ast
U
%}

A1    = 0.165;
A2    = 0.335;
b1    = 0.0455;
b2    = 0.3;
twoUC = 2*U/ch;
A32   = -b1*b2*(twoUC^2);
A33   = -(b1 + b2)*twoUC;
K1    = (A1 + A2)*b1*b2*(twoUC^2);
K2    = (A1*b1 + A2*b2)*twoUC;
K3    = 1 - A1 - A2;

tau   = 4.3*ch/U;

%{
% Special logic: tau is decreased at high angles-of-attack such that
% the dynamic and quasi-steady aoa's match over the relevant frequency
% band.  [This causes problems for the matching linearization, and it's
% probably not worth the trouble for the sorts of cases I'll be looking
% at.]
aoastp = ast(1);
aoastn = ast(2);
if (real(aq) > aoastp)
   tau = minc(maxc((1 - ((aq - aoastp)/0.1)*0.8), 0.2), 1)*tau;
elseif (real(aq) < aoastn)
   tau = minc(maxc((1 + ((aq - aoastn)/0.1)*0.8), 0.2), 1)*tau;
end
%}

taum1 = 1/tau;

%{
'------'
4.3*ch/U
tau
taum1
K1
K2
K3
%}

Amat = [-taum1 taum1*K1 taum1*K2;0 0 1;0 A32 A33];
Bmat = [taum1*K3;0;1];
dasdt = Amat*as + Bmat*aq;

dAdU = zeros(3,3);
dAdU(1,1) = -1/(4.3*ch);
dAdU(1,2) = K1/(4.3*ch) + 4*(A1 + A2)*b1*b2*twoUC/(tau*ch);
dAdU(1,3) = K2/(4.3*ch) + 2*(A1*b1 + A2*b2)/(tau*ch);
dAdU(3,2) = -8*b1*b2*U/(ch^2);
dAdU(3,3) = -2*(b1 + b2)/ch;

dBdU = zeros(3,1);
dBdU(1) = K3/(4.3*ch);

A = Amat;
By = zeros(3,2);
By(:,1) = Bmat;
By(:,2) = dAdU*as + dBdU*aq;

%'---end aoadyn---'
