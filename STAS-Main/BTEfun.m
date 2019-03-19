function [dxdt,A] = BTEfun (x,yin,params);

[dxdt,yout,A,B,C,D] = buildTurbineElectric (1,x,yin,params);

