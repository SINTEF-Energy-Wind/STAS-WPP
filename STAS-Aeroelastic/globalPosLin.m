function [xng,Dy] = globalPosLin (qB,PB,qn,Pn)
%
% Position of nodes in global coordinates.
%
%   States:           y vector:         u vector:
%                     qB
%                     qn
%
% Version:        Changes:
% --------        -------------
% 25.01.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 25.01.2018      xng matches globalPosition.m.  Derivatives verified
%                 using complex step, xng output.
%
% Inputs:
% -------
% qB,PB           : Reference nodes.
% qn,Pn           : Nodes wrt reference.
%
% Outputs:
% --------
% xng             : Vector of 3 positions.
% Dy              : State-space matrix.  3-by-12.  1:6, qB.  7:12, qn.

Dy = zeros (3,12);

TB0g = TFromTheta (PB(4:6));
[TBB0,dTBB0] = dTdth (qB(4:6));

Dy(:,1:3) = eye(3);

for jj = 1:3
   jc3 = 3*(jj-1);
   Dy(:,jj+3) = TB0g*dTBB0(:,jc3+[1:3])*(Pn(1:3) + qn(1:3));
end

TT = TB0g*TBB0;

Dy(:,7:9) = TT;
xng = PB(1:3) + qB(1:3) + TT*(Pn(1:3) + qn(1:3));
