function xng = globalPosition (qB,PB,qn,Pn)
%
% Position of nodes in global coordinates.
%
% Version:        Changes:
% --------        -------------
% 14.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 14.12.2017      Checked output for a wind turbine model.
%
% Inputs:
% -------
% qB,PB           : Reference nodes.
% qn,Pn           : Nodes wrt reference.
%
% Outputs:
% --------
% xng             : 3*Nnod vector of positions.

Nnod = size(qn,1)/6;

xng = zeros(3*Nnod,1);

for inod = 1:Nnod

   i6 = 6*(inod-1);
   i3 = 3*(inod-1);

   TB0g = TFromTheta (PB(i6+[4:6]));
   TBB0 = TFromTheta (qB(i6+[4:6]));

   xng(i3+[1:3]) = PB(i6+[1:3]) + qB(i6+[1:3]) ...
                 + TB0g*TBB0*(Pn(i6+[1:3]) + qn(i6+[1:3]));

end