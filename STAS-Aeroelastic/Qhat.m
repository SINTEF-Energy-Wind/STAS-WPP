function Qh = Qhat (qq,PP,idofs,inods)
%
% Qhat multiplies the nodal force vector to give the applied forces,
% associated with the generalized coordinates of the unconstrained
% structure.
%                        --                                                           --
%                        |                                dTB^B0              dTB^B0   |
% Qhat^T = T_B0^B T_g^B0 | (Ah*T_B0^g*T_B^B0)i + A T_B0^g ------ q + B T_B0^g ------ P |
%                        |                                  dqi                 dqi    |
%                        --                                                           --
% Version:        Changes:
% --------        -------------
% 24.11.2017      Original code.  Complete rewrite of the old version,
%                 employing calls to Qunod.m rather than buildAB.m.
%
% Version:        Verification:
% --------        -------------
% 24.11.2017      
%
% Inputs:
% -------
% qq              : Nodal DOFs.  Ref node: body relative to global.  
%                   Other nodes: node relative to undeformed, body coords.
% PP              : Undeformed nodal positions relative to body origin.
% idofs,inods     : Reference DOFs and nodes.
%

Nbod = 7;
Ndj = size(qq,1);

Qh = spalloc(Ndj,Ndj,12*Ndj);

for ibod = 1:Nbod

   if (ibod == 1)
      Nnod = inods(2) - inods(1);
      idref = idofs(1);
   elseif (ibod == 2)
      Nnod = inods(3) - inods(2);
      idref = idofs(2);
   elseif (ibod == 3)
      Nnod = inods(4) - inods(3);
      idref = idofs(3);
   elseif (ibod == 4)
      Nnod = inods(6) - inods(4);
      idref = idofs(4);
   elseif (ibod == 5)
      Nnod = inods(7) - inods(6);
      idref = idofs(6);
   elseif (ibod == 6)
      Nnod = inods(8) - inods(7);
      idref = idofs(7);
   elseif (ibod == 7)
      Nnod = inods(8) - inods(7);
      idref = idofs(8);
   end

   for inod = 1:Nnod

      noddof = idref + 6*(inod-1);

      qB = qq(idref+[1:6]);
      PB = PP(idref+[1:6]);

      if (noddof ~= idref)
         qn = qq(noddof+[1:6]);
         Pn = PP(noddof+[1:6]);
      else
         qn = zeros(6,1);
         Pn = zeros(6,1);
         Pn(4:6) = PP(noddof+6+[4:6]);
      end

      Q = Qunod (qn,qB,Pn,PB);
      Qh(idref+[1:6],noddof+[1:6]) = Q(:,1:6).';
      if (noddof ~= idref)
         Qh(noddof+[1:6],noddof+[1:6]) = Q(:,7:12).';
      end

   end


end

