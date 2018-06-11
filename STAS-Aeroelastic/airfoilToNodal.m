function F = airfoilToNodal (s,a,Fa,Tas,qq,Pn,idofs)
%
% Convert an airfoil force vector to nodal forces in pitch coordinates.
%
% Version:        Changes:
% --------        -------------
% 17.04.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 17.04.2018      
%
% Inputs:
% -------
% 
%
% Outputs:
% --------
% 

Ndj = size(qq,1);
F   = zeros(Ndj,1);

for ibod = 5:7

   if (ibod == 5)
      idref = idofs(6);
      Nel   = s.blade(1).Nel;
      conns = s.blade(1).conn;
   elseif (ibod == 6)
      idref = idofs(7);
      Nel   = s.blade(2).Nel;
      conns = s.blade(2).conn;
   elseif (ibod == 7)
      idref = idofs(8);
      Nel   = s.blade(3).Nel;
      conns = s.blade(3).conn;
   end

   for iel = 1:Nel

      ic3 = 3*Nel*(ibod-5) + 3*(iel-1);
      ic6 = 6*Nel*(ibod-5) + 6*(iel-1);

      conn  = conns(:,iel);
      rdof  = idref + 6*(conn(1)-1);
      n1dof = idref + 6*(conn(2)-1);
      n2dof = idref + 6*(conn(3)-1);

      qB  = qq(rdof+[1:6]);
      qn1 = qq(n1dof+[1:6]);
      qn2 = qq(n2dof+[1:6]);
      PB  = Pn(rdof+[1:6]);
      Pn1 = Pn(n1dof+[1:6]);
      Pn2 = Pn(n2dof+[1:6]);

      if (n1dof == rdof)

         qn1 = zeros(6,1);
         Pn1(1:3) = zeros(3,1);

         % For purposes of computing Qu's, give the reference node the
         % same undeformed orientation as the second node.
         Pn1(4:6) = Pn(n1dof+6+[4:6]);  

      end

      if (n2dof == rdof)
         qn2 = zeros(6,1);
         Pn2(1:3) = zeros(3,1);
         Pn2(4:6) = Pn(n2dof+6+[4:6]);
      end

      [xe,TsB] = elementCSFromNodes (qn1,qn2,Pn1,Pn2);
      dTsB     = derivElementCS (qn1,qn2,Pn1,Pn2,TsB);

      [Fpe,Dy] = airfoilToPitch (Fa(ic6+[1:6]),TsB,dTsB,Tas(:,ic3+[1:3]));

      F(n1dof+[1:6]) = F(n1dof+[1:6]) + 0.5*Fpe;
      F(n2dof+[1:6]) = F(n2dof+[1:6]) + 0.5*Fpe;

   end
end
