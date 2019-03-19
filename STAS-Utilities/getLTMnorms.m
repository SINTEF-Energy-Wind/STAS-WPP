function [length,time,mass,current,voltage,        ...
          velocity,force,power,stress,ndens,nvisc, ...
          stiffness,damping,resistance,inductance, ...
          capacitance,flux] = getLTMnorms (fname)

load(fname);
length      = LTMnorms(1);
time        = LTMnorms(2);
mass        = LTMnorms(3);
current     = LTMnorms(4);
velocity    = length/time;
force       = mass*length/(time^2);
power       = force*velocity;
stress      = force/(length^2);
ndens       = mass/(length^3);
nvisc       = mass/(length*time);
stiffness   = force/length;
damping     = force*time/length;
voltage     = power/current;
resistance  = voltage/current;
inductance  = voltage*time/current;
capacitance = current*time/voltage;
flux        = voltage*time;