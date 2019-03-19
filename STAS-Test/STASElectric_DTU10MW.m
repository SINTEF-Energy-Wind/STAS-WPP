function params = STASElectric_DTU10MW ()
%
% Version:        Changes:
% --------        -------------
% 03.12.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 03.12.2018      
%
% params          : --------------------- PMSG ------------------------
%                   1:    np    (-)     Number of poles.
%                   2-5:  Ig    (A)     RMS phase currents for nonlinear inductance.
%                   6-9:  Lg    (H)     Nonlinear phase inductances.
%                   10:   Rg    (Ohms)  Phase resistance.
%                   11:   lamr  (Wb)    d-component flux linkage per phase.
%                   ------- 11 --------- DC link ----------------------
%                   1:    Cdc   (F)     DC link capacitance.
%                   2:    etac  (-)     Converter efficiency.
%                   ------- 13 ------- Transformer --------------------
%                   1:    a     (-)     Np/Ns turns ratio.
%                   2-5:  It    (A)     RMS phase currents for nonlinear inductance.
%                   6-9:  Lt    (H)     Lp + (a^2)Ls nonlinear phase inductances.
%                   10:   Rt    (Ohms)  Phase resistance.
%                   ------- 23 --- Gen. current control ---------------
%                   1-4:  Ig    (A)     RMS phase currents for nonlinear inductance.
%                   5-8:  Lg    (H)     Nonlinear phase inductances.
%                   9:    lamr  (Wb)    d-component flux linkage per phase.
%                   10:   KP    (1/s)   Proportional gain.
%                   11:   KI    (1/s^2) Integral gain.
%                   12:   KF    (-)     Feed-forward gain on flux linkage.
%                   13:   ag    (rad/s) Generator current filter frequency.
%                   14:   aw    (rad/s) Generator electric speed filter frequency.
%                   ------- 37 ------------ PLL -----------------------
%                   1:    av    (rad/s) Filter on voltage measurement.
%                   2:    KPe   (rad/s) Proportional gain.
%                   3:    KIe   (rad/s2) Integral gain.
%                   4:    weh   (rad/s) Reference frequency.
%                   ------- 41 -- Reactive power control --------------
%                   1:    ap    (rad/s) ip measurement filter.
%                   2:    ais   (rad/s) is measurement filter.
%                   3:    avs   (rad/s) vs measurement filter.
%                   4:    adc   (rad/s) Vdc measurement filter.
%                   5:    KpDC  (A/V)   Proportional gain on ipd control.
%                   6:    KiDC  (A/Vs)  Integral gain on ipd control.
%                   7:    KpQ   (1/V)   Proportional gain on ipq control.
%                   8:    KiQ   (1/Vs)  Integral gain on ipq control.
%                   9:    Kppd  (1/s)   Proportional gain on vpd control.
%                  10:    Kipd  (1/s^2) Integral gain on vpd control.
%                  11:    Kppq  (1/s)   Proportional gain on vpq control.
%                  12:    Kipq  (1/s^2) Integral gain on vpq control.
%                  13:    KFp   (-)     Feedforward gain on vms in vp control.
%                  14-17: It    (A)     RMS phase currents for nonlinear inductance.
%                  18-21: Lt    (H)     Lp + (a^2)Ls nonlinear phase inductances.
%                   ------- 62 ------- Miscellaneous ------------------
%                   1:    Pnll  (W)     No-load losses.

params = zeros(63,1);

load('-ascii','LTMnorms.txt');
length = LTMnorms(1);
time   = LTMnorms(2);
mass   = LTMnorms(3);
force  = mass*length/(time^2);
power  = mass*(length^2)/(time^3);
voltage = sqrt(power);
current = power/voltage;
resistance = voltage/current;
inductance = voltage*time/current;
capacitance = current*time/voltage;
flux   = voltage*time;

twopi = 2*pi;

% ---------------------------------- PMSG --------------------------------------
np    = 198;                   % np    (-)     Number of poles.
Ig0   = 2968.4;                % (From gendes.m)
Lg0   = 0.00529;               % (From gendes.m)
Ig    = [0 1 2 3]*Ig0;         % Ig    (A)     RMS phase currents for nonlinear inductance.
Lg    = [1 1 1 1]*Lg0;         % Lg    (H)     Nonlinear phase inductances.
Rg    = 0.0366;                % Rg    (Ohms)  Phase resistance.
lamr  = 25.906;                % lamr  (Wb)    d-component flux linkage per phase.
% --------------------------------- DC Link ------------------------------------
Cdc   = 0.00437;               % Cdc   (F)     DC link capacitance.
etac  = 0.9925;                % etac  (-)     Converter efficiency.
% ------------------------------- Transformer ----------------------------------
a     = 1/8.75;                % a     (-)     Np/Ns turns ratio.
It0   = Ig0;                   % (Only a reference value for scheduling inductances.)
Lt0   = 0.0001528;             % (From gendes.m)
It    = [0 1 2 3]*It0;         % It    (A)     RMS phase currents for nonlinear inductance.
Lt    = [1 1 1 1]*Lt0;         % Lt    (H)     Lp + (a^2)Ls nonlinear phase inductances.
Rt    = 0.003;                 % Rt    (Ohms)  Phase resistance.
% --------------------------- Gen. current control -----------------------------
Ige   = Ig;                    % Ig    (A)     Estimated RMS phase currents for nonlinear inductance.
Lge   = Lg;                    % Lg    (H)     Estimated nonlinear phase inductances.
lamre = lamr;                  % lamr  (Wb)    Estimated d-component flux linkage per phase.
KP    = 50;                    % KP    (1/s)   Proportional gain.
KI    = 400;                   % KI    (1/s^2) Integral gain.
KF    = 0;                     % KF    (-)     Feed-forward gain on flux linkage.
ag    = 20*twopi;              % ag    (rad/s) Generator current filter frequency.
aw    = 20*twopi;              % aw    (rad/s) Generator electric speed filter frequency.
% ----------------------------------- PLL --------------------------------------
av    = 20*twopi;              % av    (rad/s) Filter on voltage measurement.
KPe   = 100;                   % KPe   (rad/s) Proportional gain.
KIe   = 1;                     % KIe   (rad/s2) Integral gain (should be zero, but then A is singular).
weh   = 50*twopi;              % weh   (rad/s) Reference network frequency.
% -------------------------- Reactive power control ----------------------------
ap    = 20*twopi;              % ap    (rad/s) ip measurement filter.
ais   = 20*twopi;              % ais   (rad/s) is measurement filter.
avs   = 20*twopi;              % avs   (rad/s) vs measurement filter.
adc   = 20*twopi;              % adc   (rad/s) Vdc measurement filter.
KpDC  = -0.5;                  % KpDC  (A/V)   Proportional gain on ipd control.
KiDC  = -0.1;                  % KiDC  (A/Vs)  Integral gain on ipd control.
KpQ   = -0.00002;              % KpQ   (1/V)   Proportional gain on ipq control.
KiQ   = -0.006;                % KiQ   (1/Vs)  Integral gain on ipq control.
Kppd  = 600;                   % Kppd  (1/s)   Proportional gain on vpd control.
Kipd  = 1000;                  % Kipd  (1/s^2) Integral gain on vpd control.
Kppq  = 600;                   % Kppq  (1/s)   Proportional gain on vpq control.
Kipq  = 1000;                  % Kipq  (1/s^2) Integral gain on vpq control.
KFp   = 0;                     % KFp   (-)     Feedforward gain on vms in vp control.
Ite   = It;                    % It    (A)     Estimated RMS phase currents for nonlinear inductance.
Lte   = Lt;                    % Lt    (H)     Estimated Lp + (a^2)Ls nonlinear phase inductances.
% ------------------------------ Miscellaneous --------------------------------
Pnnl  = 7.33e4;                % Pnll  (W)     No-load (constant) losses.

% ---------------------------------- PMSG --------------------------------------
params(1)     = np;
params(2:5)   = Ig      / current;
params(6:9)   = Lg      / inductance;
params(10)    = Rg      / resistance;
params(11)    = lamr    / flux;
% --------------------------------- DC Link ------------------------------------
params(12)    = Cdc     / capacitance;
params(13)    = etac;       
% ------------------------------- Transformer ----------------------------------
params(14)    = a;
params(15:18) = It      / current;
params(19:22) = Lt      / inductance;
params(23)    = Rt      / resistance;
% --------------------------- Gen. current control -----------------------------
params(24:27) = Ige     / current;
params(28:31) = Lge     / inductance;
params(32)    = lamre   / flux;
params(33)    = KP      * time;
params(34)    = KI      * time^2;
params(35)    = KF;
params(36)    = ag      * time;
params(37)    = aw      * time;
% ----------------------------------- PLL --------------------------------------
params(38)    = av      * time;
params(39)    = KPe     * time;
params(40)    = KIe     * time^2;
params(41)    = weh     * time;
% -------------------------- Reactive power control ----------------------------
params(42)    = ap      * time;
params(43)    = ais     * time;
params(44)    = avs     * time;
params(45)    = adc     * time;
params(46)    = KpDC    / (current/voltage);
params(47)    = KiDC    / (current/(voltage*time));
params(48)    = KpQ     * voltage;
params(49)    = KiQ     * voltage*time;
params(50)    = Kppd    * time;
params(51)    = Kipd    * time^2;
params(52)    = Kppq    * time;
params(53)    = Kipq    * time^2;
params(54)    = KFp;
params(55:58) = Ite     / current;
params(59:62) = Lte     / inductance;
% ------------------------------ Miscellaneous ---------------------------------
params(63)    = Pnnl    / power;