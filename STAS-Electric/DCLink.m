function [dxdt,yout,A,By,C,Dy] = DCLink (Linflag,x,yin,params)
%
% An equivalent-circuit model of the converters and DC link, with
% converters represented only by an efficiency.  The governing 
% equation of the DC link is
% dVdc/dt = (1/Cdc) (IIn - IIg)
%
%   states:           y vector:
%   Vdc        1      igd,q     1,2    in   (generator)
%                     vgd,q     3,4    in   (converter control)
%                     ipd,q     5,6    in   (transformer)
%                     vpd,q     7,8    in   (converter control)
%                     IIg        9     out
%                     IIn       10     out
%
% Vdc is the DC link voltage.  IIg is the current fed to the DC link
% by the generator-side converter.  IIn is the current extracted by
% the network-side converter.  ig is the AC generator current in the
% dq frame, vg is the AC generator voltage, ip is the AC current at
% the transformer primary terminals, vp is the AC voltage at the
% transformer primary terminals.
%
% Version:        Changes:
% --------        -------------
% 02.05.2016      Original code.
% 11.07.2017      Adapted for complex step derivatives.
% 27.08.2018      Modified for linear/nonlinear equation pairs.
%
% Version:        Verification:
% --------        -------------
% 02.05.2016      
% 11.07.2017      
% 27.08.2018      
%
% Inputs:
% -------
% Linflag         : set to 1 to perform linearization.
% x               : 1:   Vdc   (V)   DC link voltage
% yin             : 1,2: ig    (A)   stator current, dq frame
%                   3,4: vg    (V)   stator voltage, dq frame
%                   5,6: ip    (A)   transformer-side current, dq frame
%                   7,8: vp    (V)   transformer-side voltage, dq frame
% params          : 1:   Cdc   (F)   DC link capacitance.
%                   2:   etac  (-)   Converter efficiency.
%
% Outputs:
% --------
% A, By, C, Dy

Nx = 1;
Ny = 10;

dxdt = 0;
yout = zeros(2,1);

A  = sparse(Nx,Nx);
By = sparse(Nx,Ny);
C  = sparse(Ny,Nx);
Dy = sparse(Ny,Ny);

Vdc = x;

ig = yin(1:2);
vg = yin(3:4);
ip = yin(5:6);
vp = yin(7:8);

Cdc = params(1);
eta = params(2);

IIg = eta*(ig.')*vg/Vdc;
IIn = (ip.')*vp/(eta*Vdc);

yout(1) = IIg;
yout(2) = IIn;
dxdt = (IIg - IIn)/Cdc;

if (Linflag == 1)

   % DC link dynamics.
   L = 1;
   By(9)  = 1/Cdc;
   By(10) = -By(9);

   % IIg current equation.
   C(9) = -IIg/Vdc;
   Dy(9,1:2) = eta*(vg.')/Vdc;
   Dy(9,3:4) = eta*(ig.')/Vdc;

   % IIn current equation.
   C(10) = -IIn/Vdc;
   Dy(10,5:6) = (vp.')/(eta*Vdc);
   Dy(10,7:8) = (ip.')/(eta*Vdc);

end

