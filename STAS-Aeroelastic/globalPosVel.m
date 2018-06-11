function [xng,vng,Dy] = globalPosVel (s,idofs,qq,PP,dqqdt)
%
% Linearized global positions and velocities for all the nodes in the
% structure.
%
%   States:           y vector:         u vector:
%                     q
%                     dq/dt
%
% Version:        Changes:
% --------        -------------
% 25.01.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 25.01.2018      Checked complex step gradients using xng,vng outputs.
%
% Inputs:
% -------
% qq,dqq/dt       : Unconstrained, un-reduced nodal DOFs.
% PP              : Undeformed pos,rot of nodes.
%
% Outputs:
% --------
% xng,vng         : Nodal positions and velocities at the input conditions.
% Dy              : State matrix.

Ndof = size(qq,1);
Nn = Ndof/6;

Dy = spalloc (3*Nn+Ndof,2*Ndof,round(0.05*Ndof*Ndof));
xng = zeros(3*Nn,1);
vng = zeros(Ndof,1);

ix = 0;
for ibod = 1:7

   if (ibod == 1)
      idref = idofs(1);
      Nnod  = s.foundation.Nnod;
   elseif (ibod == 2)
      idref = idofs(2);
      Nnod  = s.tower.Nnod;
   elseif (ibod == 3)
      idref = idofs(3);
      Nnod  = s.nacelle.Nnod;
   elseif (ibod == 4)
      idref = idofs(4);
      Nnod   = s.driveshaft.Nnod;
   elseif (ibod == 5)
      idref = idofs(6);
      Nnod   = s.blade(1).Nnod;
   elseif (ibod == 6)
      idref = idofs(7);
      Nnod   = s.blade(2).Nnod;
   elseif (ibod == 7)
      idref = idofs(8);
      Nnod   = s.blade(3).Nnod;
   end

   for inod = 1:Nnod

      rdof = idref;
      mdof = idref + 6*(inod-1);

      qB    = qq(rdof+[1:6]);
      PB    = PP(rdof+[1:6]);
      dqBdt = dqqdt(rdof+[1:6]);

      if (mdof == rdof)
         qn    = zeros(6,1);
         Pn    = zeros(6,1);
         dqndt = zeros(6,1);
         Pn(4:6) = PP(rdof+6+[4:6]);  % Consistent with other functions.
      else
         qn    = qq(mdof+[1:6]);
         Pn    = PP(mdof+[1:6]);
         dqndt = dqqdt(mdof+[1:6]);
      end

      [xn,ddy] = globalPosLin (qB,PB,qn,Pn);
      xng(ix+[1:3]) = xn;
      Dy(ix+[1:3],rdof+[1:6]) = ddy(:,1:6);
      if (mdof ~= rdof)
         Dy(ix+[1:3],mdof+[1:6]) = ddy(:,7:12);
      end

      [vn,ddy] = globalVelocity (qn,qB,Pn,PB,dqndt,dqBdt);
      vng(mdof+[1:6]) = vn;
      Dy(3*Nn+mdof+[1:6],rdof+[1:6])      = ddy(:,1:6);
      Dy(3*Nn+mdof+[1:6],Ndof+rdof+[1:6]) = ddy(:,13:18);

      if (mdof ~= rdof)
         Dy(3*Nn+mdof+[1:6],mdof+[1:6])      = ddy(:,7:12);
         Dy(3*Nn+mdof+[1:6],Ndof+mdof+[1:6]) = ddy(:,19:24);
      end

      ix = ix + 3;

   end

end