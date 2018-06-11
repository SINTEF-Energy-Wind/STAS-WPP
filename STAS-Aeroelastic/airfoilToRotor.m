function [Tar,dTar] = airfoilToRotor (qy,qp,qn1,qn2,Py,Pp,Pn1,Pn2,Try,Tas)
%
% The transform from airfoil to rotor coordinates.
%
% Version:        Changes:
% --------        -------------
% 21.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 21.12.2017      dTar verified by complex step.
%
% Inputs:
% -------
% qy,Py           : Nacelle body reference.
% qp,Pp           : Blade body reference.
% qn,Pn           : Nodes.
% Try,Tas         : Constant matrices.
%
% Outputs:
% --------
% Tar             : Transform from airfoil to rotor coordinates.
% dTar            : Derivatives wrt Phiy, Phip, qn1, qn2.

dTar = zeros (3,3*24);  % 1:6: qy.  7:12: qp.  13:18: qn1.  19:24: qn2.

[Tyy0,dTyy0] = dTdth (qy(4:6));
[Tpp0,dTpp0] = dTdth (qp(4:6));
[xe,Tsp]     = elementCSFromNodes (qn1,qn2,Pn1,Pn2);
dTsp         = derivElementCS (qn1,qn2,Pn1,Pn2,Tsp);
Ty0g         = TFromTheta (Py(4:6));
Tp0g         = TFromTheta (Pp(4:6));
Tag          = Tp0g*Tpp0*Tsp*Tas;
Tgr          = (Ty0g*Tyy0*Try).';

Tar = Tgr*Tag;

for jj = 1:3
   jc3a = 3*(jj+3-1);
   jc3b = 3*(jj-1);
   dTar(:,jc3a+[1:3]) = ((Ty0g*dTyy0(:,jc3b+[1:3])*Try).')*Tag;
end
for jj = 1:3
   jc3a = 3*(jj+9-1);
   jc3b = 3*(jj-1);
   dTar(:,jc3a+[1:3]) = Tgr*Tp0g*dTpp0(:,jc3b+[1:3])*Tsp*Tas;
end
for jj = 1:12
   jc3a = 3*(jj+12-1);
   jc3b = 3*(jj-1);
   dTar(:,jc3a+[1:3]) = Tgr*Tp0g*Tpp0*dTsp(:,jc3b+[1:3])*Tas;
end