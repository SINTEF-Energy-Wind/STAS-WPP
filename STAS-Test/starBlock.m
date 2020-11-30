function [eta,D] = starBlock (u,uLB,uUB,eps,weta,linFlag)
%
% Soft saturation block.
%
% Version:        Changes:
% --------        -------------
% 19.10.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 19.10.2020      Sample output checked.  Derivatives verified using
%                 complex step.
%
% Inputs:
% -------
% u,uLB,uUB       : The saturating signal and its saturation bounds.
% eps             : The error input to the PI.
% weta            : Width of the transition region in the saturations.
% linFlag         : Set to 1 if linearized D output is desired.
%
% Outputs:
% --------
% eta             : Anti-windup output.
% D               : = d eta / d[u, uLB, uUB, eps].


del = minc (uUB - u, u - uLB);          % minc: for complex step.
kap = sign (real(u - (uUB + uLB)/2));   % real: for complex step.

p = [0;1 + 0.1/weta];
[ye,dye,d2ye] = saturate (-(kap*eps)/weta,p);
[yd,dyd,d2yd] = saturate (del/weta,p);

ne = 0.5*(1 + ye);
nd = 0.5*(1 + yd);

eta = nd + ne*(1 - nd);

%[kap eps -(kap*eps)/weta ye del del/weta yd ne nd eta]

D = zeros(1,4);
if (linFlag == 1)

   detadnd = 1 - ne;
   detadne = 1 - nd;
   dnddyd = 0.5;
   dnedye = 0.5;

   if (real(uUB-u) <= real(u-uLB))      % real: for complex step.
      ddeldu = -1;
      ddelduUB = 1;
      ddelduLB = 0;
   else
      ddeldu = 1;
      ddelduUB = 0;
      ddelduLB = -1;
   end

   ddd = detadnd*dnddyd*dyd;

   D(1) = ddd*ddeldu/weta;
   D(2) = ddd*ddelduLB/weta;
   D(3) = ddd*ddelduUB/weta;
   D(4) = -detadne*dnedye*dye*kap/weta;

end
