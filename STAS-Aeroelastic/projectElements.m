function [y,Dia,Dy,Dd] = projectElements (q,P,iq,idofs,idofm,Tn_y)
%
% Construct the relationships between projected rotorplane
% coordinates and q DOFs.
%
% Output y's (rows of Dy) are xeg, xn1, xn2, xer, r, and Lp.
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
%              (Repeat for each element.)
%                     Dp     (3)
%
% Version:        Changes:
% --------        -------------
% 08.02.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 08.02.2018      Partial verification with complex step via qy,qB,qn1,qn2,xhg
%                 inputs.
%
% Inputs:
% -------
% iq              : Indices from BEMprep.m.
%
% Outputs:
% --------
% 

Nel = size (iq,2);
Neb = Nel/3;

Dy = spalloc (14*Nel,45*Nel+3,60*Nel);
Dd = spalloc (3,9,9);  % 9: xnr at tip nodes.
y  = zeros(14*Nel,1);

Try = Tn_y;
Py = P(iq(1:6,1));
qy = q(iq(1:6,1));
Ty0g = TFromTheta (Py(4:6));
[Tyy0,dTyy0] = dTdth (qy(4:6));
Tgr = (Ty0g*Tyy0*Try).';
xhg  = globalPosition (q(idofs(4)+[1:6]),P(idofs(4)+[1:6]), ...
                       q(idofm(6)-6+[1:6]),P(idofm(6)-6+[1:6]));

for iel = 1:Nel

   ic131 = 131*(iel-1);
   ic45  =  45*(iel-1);
   ic14  =  14*(iel-1);

   qy  = q(iq(1:6,iel));
   qB  = q(iq(7:12,iel));
   qn1 = q(iq(13:18,iel));
   qn2 = q(iq(19:24,iel));

   Py  = P(iq(1:6,iel));
   PB  = P(iq(7:12,iel));
   Pn1 = P(iq(13:18,iel));
   Pn2 = P(iq(19:24,iel));

   if (iq(7,iel) == iq(13,iel))  % Node 1 is the reference.
      qn1 = zeros(6,1);
      Pn1(1:3) = zeros(3,1);
      Pn1(4:6) = Pn2(4:6);
   end

   [xeg,xnr1,xnr2,xer,r,Lp,ddp] =  ...
          projectElement (qy,qB,qn1,qn2,Py,PB,Pn1,Pn2,Try,xhg);
   Dy(ic14+[1:14],ic45+[1:45]) = ddp;
   y(ic14+[1:3])   = xeg;
   y(ic14+[4:6])   = xnr1;
   y(ic14+[7:9])   = xnr2;
   y(ic14+[10:12]) = xer;
   y(ic14+[13])    = r;
   y(ic14+[14])    = Lp;

end

% Outer diameter.  Use blade-by-blade.
Dia = 0;
for ib = 1:3

   ic3 = 3*(ib-1);

   Bref = idofs(5+ib);
   tref = Bref + 6*Neb;
   xtg = globalPosition (q(Bref+[1:6]),P(Bref+[1:6]), ...
                         q(tref+[1:6]),P(tref+[1:6]));
   rg = xtg - xhg;
   rr = Tgr*rg;
   [Do,ddy] = Douter (rr);
   Dia = Dia + Do/3;
   Dd(ib,ic3+[1:3]) = ddy;

%ib
%Do
%ddy
%Dia

end

%Dd
