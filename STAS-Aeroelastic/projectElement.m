function [xeg,xnr1,xnr2,xer,r,Lp,Dy] = projectElement (qy,qB,qn1,qn2,Py,PB,Pn1,Pn2,Try,xhg)
%
% Project the blade element onto the rotor plane.
%
%   States:           y vector:         u vector:
%                     qy     1:6
%                     qB     7:12
%                     qn1   13:18
%                     qn2   19:24
%                     xng1  25:27
%                     xng2  28:30
%                     xeg   31:33
%                     xhg   34:36
%                     xnr1  37:39
%                     xnr2  40:42
%                     xer   43:45
%
% Version:        Changes:
% --------        -------------
% 26.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 26.12.2017      Partial verification with complex step via qy,qB,qn1,qn2,xhg
%                 inputs.
%
% Inputs:
% -------
% qy,Py           : Nacelle body reference.
% qp,Pp           : Blade body reference.
% qn,Pn           : Nodes.
% Try             : Rotor-to-yaw CS transform.
% xhg             : Position of the rotorplane center (= hub node) in 
%                   global coordinates.
%
% Outputs:
% --------
% xeg,xer         : Position of the element centroid in global and rotor
%                   coordinates.
% xnr1,2          : Position of nodes in rotor coordinates.
% r               : Radius of the element (distance from axis of rotation).
% Lp              : Projected length of the element.

Dy = zeros (14,45); % 1:3=xeg, 4:6=xnr1, 7:9=xnr2, 10:12=xer, 13=r, 14=Lp.

iqy   =  [1:6].';
iqB   =  [7:12].';
iqn1  = [13:18].';
iqn2  = [19:24].';
ixng1 = [25:27].';
ixng2 = [28:30].';
ixeg  = [31:33].';
ixhg  = [34:36].';
ixnr1 = [37:39].';
ixnr2 = [40:42].';
ixer  = [43:45].';

Ty0g = TFromTheta (Py(4:6));
[Tyy0,dTyy0] = dTdth (qy(4:6));

xng = globalPosition ([qB;qB],[PB;PB],[qn1;qn2],[Pn1;Pn2]);

xeg = 0.5*(xng(1:3) + xng(4:6));

Tgr = (Ty0g*Tyy0*Try).';

xer  = Tgr*(xeg - xhg);
xnr1 = Tgr*(xng(1:3) - xhg);
xnr2 = Tgr*(xng(4:6) - xhg);
dx = xnr2 - xnr1;

r  = sqrt(xer(1)^2 + xer(2)^2);
Lp = sqrt(dx(1)^2 + dx(2)^2);

Dy(1:3,ixng1)   =  0.5*speye(3);
Dy(1:3,ixng2)   =  0.5*speye(3);

Dy(4:6,ixng1)   =  Tgr;
Dy(4:6,ixhg)    = -Tgr;

Dy(7:9,ixng2)   =  Tgr;
Dy(7:9,ixhg)    = -Tgr;

Dy(10:12,ixeg)  =  Tgr;
Dy(10:12,ixhg)  = -Tgr;

for jj = 1:3
   jc3 = 3*(jj-1);
   mat = (Ty0g*dTyy0(:,jc3+[1:3])*Try).';
   Dy(4:6,iqy(jj+3))   = mat*(xng(1:3) - xhg);
   Dy(7:9,iqy(jj+3))   = mat*(xng(4:6) - xhg);
   Dy(10:12,iqy(jj+3)) = mat*(xeg      - xhg);
end
Dy(13,ixer)  =  [xer(1) xer(2) 0]/r;
Dy(14,ixnr2) =  [dx(1)   dx(2) 0]/Lp;  %  (dx.')/Lp; (wrong, project onto
Dy(14,ixnr1) = -[dx(1)   dx(2) 0]/Lp;  % -(dx.')/Lp;  rotorplane)

%{
'---projectElement---'
[qy qB qn1 qn2 Py PB Pn1 Pn2]
xhg
xng
xeg
Trg = (Tgr.')
xer
xnr1
xnr2
dx
r
Lp
'--end projectElement--'
%}
