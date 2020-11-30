function [Pow,y,Du,TBg] = appliedPower (s,qq,dqqdt,F,PP,idofs)
%
% Compute the overall power (force-times-velocity) of the applied 
% nodal forces. 
%
%   States:           y vector:         u vector:
%                     vg       Ndj      qB          Ndj
%                     Fg       Ndj      dqB/dt      Ndj
%                                       FB          Ndj
%
% Version:        Changes:
% --------        -------------
% 25.12.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 25.12.2019      
%
% Inputs:
% -------
% 
%
% Outputs:
% --------
% Pow             : Mean power.
% y               : v^g and F^g.
% Du              : = dy/du

Ndj = size(qq,1);
Nu = 3*Ndj;
Ny = 2*Ndj;

Du = sparse(Ny,Nu);
TBg = sparse(Ndj,Ndj);

vg0 = zeros(Ndj,1);
Fg0 = zeros(Ndj,1);
Dvg = sparse(Ndj,2*Ndj);   % 2*Ndj: q, dq/dt.
DFg = sparse(Ndj,2*Ndj);   % 2*Ndj: q, FB.

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

      [vg,ddy] = globalVelocity (qn,qB,Pn,PB,dqndt,dqBdt);
      vg0(mdof+[1:6]) = vg;

      if (mdof == rdof)
         Dvg(rdof+[1:6],[rdof+[1:6] Ndj+rdof+[1:6]]) = ddy(:,[[1:6] [13:18]]);
      else
         Dvg(mdof+[1:6],[rdof+[1:6] mdof+[1:6] Ndj+rdof+[1:6] Ndj+mdof+[1:6]]) = ddy;
      end

      [TBB0,dTBB0] = dTdth (qB(4:6));
      TB0g = TFromTheta (PB(4:6));
      T = TB0g*TBB0;
      dT = TB0g*dTBB0;
      TT = [T, zeros(3,3);zeros(3,3), T];
      Fg0(mdof+[1:6]) = TT*F(mdof+[1:6]);

      TBg(mdof+[1:3],mdof+[1:3]) = T;
      TBg(mdof+[4:6],mdof+[4:6]) = T;
      
      DFg(mdof+[1:6],Ndj+mdof+[1:6]) = TT;
      for jj = 1:3
         jc3 = 3*(jj-1);
         dTT = [dT(:,jc3+[1:3]), zeros(3,3);zeros(3,3), dT(:,jc3+[1:3])];
         DFg(mdof+[1:6],rdof+3+jj) = dTT*F(mdof+[1:6]);
      end

   end

end

TBg(Ndj-6+[1:6],Ndj-6+[1:6]) = speye(6);

y = [vg0;Fg0];
Du(1:Ndj,1:2*Ndj) = Dvg;
Du(Ndj+[1:Ndj],[[1:Ndj] 2*Ndj+[1:Ndj]]) = DFg;

Pow = sum(Fg0.*vg0);
