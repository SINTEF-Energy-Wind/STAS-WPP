function [dxdt,yout,A,By,C,Dy] = converterControlNet (Linflag,x,yin,params)
%
% The network-side converter is controlled so as to provide a commanded
% reactive power, with the active power controlled so as to hold the DC
% link voltage at its nominal value.
%
% Note, could implement a PLL instead of the simple measurement delay on
% the electrical frequency as implemented presently.  For this, let the
% pll own vms (?) and move wem to the y vector as input.
%
%   states:           y vector:
%   impd,q    1,2     ipd,q     1,2    in   (transformer)
%   imsd,q    3,4     isd,q     3,4    in   (transformer)
%   vmsd,q    5,6     vsd,q     5,6    in   (grid)
%   Vmdc       7      Vdc        7     in   (converter)
%   PsiDC      8      Vhdc       8     in   (input)
%   PsiQ       9      Qh         9     in   (voltage control)
%   PsiP     10,11    wem       10     in   (PLL)
%                     vpd,q    11,12   out
%
% Version:        Changes:
% --------        -------------
% 21.11.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 21.11.2018      Derivatives verified by complex step.
%
% Inputs:
% -------
% Linflag         : Set to 1 for linearized equations.
% x               : 1,2   imp   (A)     Measured converter terminal current.
%                   3,4   ims   (A)     Measured transformer terminal current.
%                   5,6   vms   (V)     Measured transformer terminal voltage.
%                    7    Vmdc  (V)     Measured DC link voltage.
%                    8    PsiDC (A)     Integrated DC link voltage error (w/gain).
%                    9    PsiQ  (A)     Integrated reactive power error (w/gain).
%                  10,11  PsiP  (A/s)   Integrated current error (w/gain).
% yin             : 1,2   ip    (A)     Converter terminal current.
%                   3,4   is    (A)     Transformer terminal current.
%                   5,6   vs    (A)     Transformer terminal voltage.
%                    7    Vdc   (V)     DC link voltage.
%                    8    Vhdc  (V)     DC link voltage command.
%                    9    Qh    (VA)    Reactive power command.
%                   10    wem   (rad/s) Measured electrical frequency.
% params          :  1    ap    (rad/s) ip measurement filter.
%                    2    ais   (rad/s) is measurement filter.
%                    3    avs   (rad/s) vs measurement filter.
%                    4    adc   (rad/s) Vdc measurement filter.
%                    5    KpDC  (A/V)   Proportional gain on ipd control.
%                    6    KiDC  (A/Vs)  Integral gain on ipd control.
%                    7    KpQ   (1/V)   Proportional gain on ipq control.
%                    8    KiQ   (1/Vs)  Integral gain on ipq control.
%                    9    Kppd  (1/s)   Proportional gain on vpd control.
%                   10    Kipd  (1/s^2) Integral gain on vpd control.
%                   11    Kppq  (1/s)   Proportional gain on vpq control.
%                   12    Kipq  (1/s^2) Integral gain on vpq control.
%                   13    KFp   (-)     Feedforward gain on vms in vp control.
%                  14-17  It    (A)     RMS phase currents for nonlinear inductance
%                  18-21  Lt    (H)     Lp + (a^2)Ls nonlinear phase inductances
%
% Outputs:
% --------
% yout            : 1,2   vp    (V)     Converter terminal voltage.

Nx    = size(x,1);
Nyi   = size(yin,1);
Nyo   = 2;
Ny    = Nyi + Nyo;

dxdt  = zeros(Nx,1);
yout  = zeros(Nyo,1);
A     = zeros(Nx,Nx);
By    = zeros(Nx,Ny);
C     = zeros(Ny,Nx);
Dy    = zeros(Ny,Ny);

sq3 = 1/sqrt(3);

imp   = x(1:2);
ims   = x(3:4);
vms   = x(5:6);
Vmdc  = x(7);
PsiDC = x(8);
PsiQ  = x(9);
PsiP  = x(10:11);
ip    = yin(1:2);
is    = yin(3:4);
vs    = yin(5:6);
Vdc   = yin(7);
Vhdc  = yin(8);
Qh    = yin(9);
wem   = yin(10);
ap    = params(1);
ais   = params(2);
avs   = params(3);
adc   = params(4);
KpDC  = params(5);
KiDC  = params(6);
KpQ   = params(7);
KiQ   = params(8);
Kppd  = params(9);
Kipd  = params(10);
Kppq  = params(11);
Kipq  = params(12);
KFp   = params(13);
IIt   = params(14:17);
LLt   = params(18:21);

Qms = vms(2)*ims(1) - vms(1)*ims(2);
edc = Vhdc - Vmdc;
eQ  = Qh - Qms;
ipd = KpDC*edc + PsiDC;
ipq = KpQ*eQ + PsiQ;
II  = maxc(sq3*sqrt(imp(1)^2 + imp(2)^2),eps);
[L,dL] = inductance (IIt,LLt,II);

dxdt(1:2)   = -ap*imp + ap*ip;
dxdt(3:4)   = -ais*ims + ais*is;
dxdt(5:6)   = -avs*vms + avs*vs;
dxdt(7)     = -adc*Vmdc + adc*Vdc;
dxdt(8)     =  KiDC*edc;
dxdt(9)     =  KiQ*eQ;
dxdt(10:11) =  [Kipd 0;0 Kipq]*[ipd-imp(1);ipq-imp(2)];

term = [Kppd 0;0 Kppq]*[ipd-imp(1);ipq-imp(2)] + PsiP + wem*[0 -1;1 0]*imp + KFp*vms;
yout = L*term;

if (Linflag == 1)

   A(1:2,1:2)    =  [-ap 0;0 -ap];
   A(3:4,3:4)    =  [-ais 0;0 -ais];
   A(5:6,5:6)    =  [-avs 0;0 -avs];
   A(7,7)        = -adc;
   A(8,7)        = -KiDC;
   A(9,3:6)      = -KiQ*[vms(2), -vms(1), -ims(2), ims(1)];
   A(10,1)       = -Kipd;
   A(10,7)       = -Kipd*KpDC;
   A(10,8)       =  Kipd;
   A(11,2)       = -Kipq;
   A(11,3:6)     = -Kipq*KpQ*[vms(2), -vms(1), -ims(2), ims(1)];
   A(11,9)       =  Kipq;

   By(1:2,1:2)   =  [ap 0;0 ap];
   By(3:4,3:4)   =  [ais 0;0 ais];
   By(5:6,5:6)   =  [avs 0;0 avs];
   By(7,7)       =  adc;
   By(8,8)       =  KiDC;
   By(9,9)       =  KiQ;
   By(10,8)      =  Kipd*KpDC;
   By(11,9)      =  Kipq*KpQ;

   C(Nyi+[1:2],1:2) =  L*wem*[0 -1;1 0]                 ...
                    +  dL*term*[imp(1)/II, imp(2)/II]/3 ...
                    -  L*[Kppd, 0;0, Kppq];
   C(Nyi+[1:2],5:6) =  L*[KFp 0;0 KFp];
   C(Nyi+1,7)    = -L*Kppd*KpDC;
   C(Nyi+1,8)    =  L*Kppd;
   C(Nyi+2,3:6)  =  C(Nyi+2,3:6) - L*Kppq*KpQ*[vms(2), -vms(1), -ims(2), ims(1)];
   C(Nyi+2,9)    =  L*Kppq;
   C(Nyi+[1:2],10:11) = L*eye(2);

   Dy(Nyi+1,8)   =  L*Kppd*KpDC;
   Dy(Nyi+2,9)   =  L*Kppq*KpQ;
   Dy(Nyi+[1:2],10) =  L*[0 -1;1 0]*imp;

end

