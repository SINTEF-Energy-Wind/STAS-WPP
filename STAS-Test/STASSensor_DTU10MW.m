function m = STASSensor_DTU10MW ()
%
% Version:        Changes:
% --------        -------------
% 15.02.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 15.02.2020      
%
% params          : 1: aW    (rad/s)    W    time constant.
%                   2: ab    (rad/s)    beta time constant.
%                   3: ay    (rad/s)    yaw  time constant.
%                   4: aP    (rad/s)    Pe   time constant.
%                   5: ad    (rad/s)    dnac time constant.
%                   6: aV    (rad/s)    V    time constant.
%                   7: at    (rad/s)    thw  time constant.

load 'LTMnorms.txt';
[length,time,mass,current,voltage,        ...
 velocity,force,power,stress,ndens,nvisc, ...
 stiffness,damping,resistance,inductance, ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

twopi = 2*pi;

% A bit slow, can speed them up if this causes problems.
m.aW = 2     *twopi    *time;
m.ab = 2     *twopi    *time;
m.ay = 0.5   *twopi    *time;
m.aP = 4     *twopi    *time;
m.ad = 2     *twopi    *time;
m.aV = 0.10  *twopi    *time;
m.at = 0.10  *twopi    *time;
