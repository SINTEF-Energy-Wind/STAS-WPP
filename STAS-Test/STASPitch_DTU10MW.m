function params = STASPitch_DTU10MW ()
%
% Version:        Changes:
% --------        -------------
% 21.01.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 21.01.2019      
%
% params          : 1: ka    (Nm/rad)   Torsional stiffness.
%                   2: ca    (Nms/rad)  Torsional damping.
%                   3: aa    (rad/s)    Frequency of motor dynamics.
%                   4: K     (1/s)      Gain on pitch rate control.
%                   5: bmax  (rad)      Max slew angle.
%                   6: bmin  (rad)      Min slew angle.
%                   7: bdmax (rad/s)    Max slew rate.

load 'LTMnorms.txt';
[length,time,mass,current,voltage,        ...
 velocity,force,power,stress,ndens,nvisc, ...
 stiffness,damping,resistance,inductance, ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

p180 = pi/180;

ka    = 1e8              / (force*length);
ca    = 5e6              / (force*length*time);
aa    = 1*(2*pi)         * time;
K     = 4                * time;
bmax  = 90*p180;
bmin  = -5*p180;
bdmax =  9*p180;

params = [ka;ca;aa;K;bmax;bmin;bdmax];
