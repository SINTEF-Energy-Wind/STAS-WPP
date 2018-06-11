function [T,dT,d2T] = bodyTransforms (idofs,qB)
%
% This computes the TB^B0 transforms based on an input q vector.
%
% Version:        Changes:
% --------        -------------
% 16.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 16.10.2017      Based on dTdth and d2Tdth2, which have been verified.
%

Nbod = 7;

T = zeros(3,3*Nbod);
dT = zeros(3,3*3*Nbod);
d2T = zeros(3,3*3*3*Nbod);

ind = [1 2 3 4 6 7 8];

for ibod = 1:7

   ic3 = 3*(ibod-1);
   ic9 = 9*(ibod-1);
   ic27 = 27*(ibod-1);

   Th = qB(idofs(ind(ibod))+[4:6]);
   [T(:,ic3+[1:3]),dT(:,ic9+[1:9])] = dTdth(Th);
   d2T(:,ic27+[1:27]) = d2Tdth2 (Th,T(:,ic3+[1:3]),dT(:,ic9+[1:9]));

end
