function [dxdt,A,B] = shaftFunc (x,u,p,cpct,KeTab)
%
% Version:        Changes:
% --------        -------------
% 08.02.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 08.02.2019      
%

W    = x(1);
Vinf = x(2);
Wm   = u(1);
bet  = u(2);
Pe   = u(3);
R    = p(1);
J    = p(2);
rho  = p(3);
KW   = p(4);
KV   = p(5);

[etag,detag] = gains1 (u(3),KeTab);

pA2  = 0.5*rho*pi*(R^2);
V3   = Vinf^3;

[cp,ct,dcp,dct] = cpvwb (cpct,R,Vinf,W,bet);

if (W < 0.1)
   dxdt = 0;
   A = zeros(2,2);
   B = zeros(2,3);
else
   term = Pe/(etag*W);
   dxdt = [cp*pA2*V3/W - term; 0]/J + [KW;KV]*(u(1) - x(1));
   A    = [pA2*V3*(dcp(2)/W-cp/(W^2))+term/W, ...
           (pA2/W)*(dcp(1)*V3+3*cp*(Vinf^2)); ...
           0, 0]/J - [KW, 0;KV, 0];
   B    = [KW, pA2*V3*dcp(3)/(W*J), -1/(etag*W*J)+detag*term/(etag*J); ...
           KV, 0, 0];
end
