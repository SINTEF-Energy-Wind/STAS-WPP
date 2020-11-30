function dQhF = dQhdqF (qq,PP,FF,idofs,inods)
%
% In the linearized equations, an effective stiffness is contributed by
% -(dQh/dq*F0).  This returns dQh/dq*F0, mark that a negative sign will
% be needed.
%
% Version:        Changes:
% --------        -------------
% 24.11.2017      Original code.
% 27.01.2020      Fixed an error, adding the "if" statement before
%                 the second idof loop.
%
% Version:        Verification:
% --------        -------------
% 24.11.2017      Linearization verified by complex step.
% 27.01.2020      Linearization verified by complex step on all DOFs.
%
% Inputs:
% -------
% qq              : Nodal DOFs.  Ref node: body relative to global.  
%                   Other nodes: node relative to undeformed, body coords.
% PP              : Undeformed nodal positions relative to body origin.
% FF              : Nodal force vector.
% idofs...inodm   : Reference DOFs and nodes.
%

Nbod = 7;
Ndj = size(qq,1);
Ndof = Ndj - 6;

dQhF = spalloc(Ndj,Ndj,12*Ndj);

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

   % For each DOF, a vector dQh/dq*F needs to be created.  This nominally
   % spans all the DOFs, but due to the separate nature of each body,
   % prior to constraints, we can consider each body in turn.  To avoid
   % computing dQu/dq multiple times for each node, we index the outer
   % loop according to the column of dQu/dq, considering 6 DOFs (one node)
   % at a time.
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
      Fn = FF(noddof+[1:6]);

      dQ = dQudq (qn,qB,Pn,PB);

      % dQ is 6-by-12*12.  The transpose is used in the present equations.
      % That is, the 6 determines which components of F are multiplied
      % (columns of Q), while the 12 corresponds to the rows of Q.  The
      % DERIVATIVE index determines the column of Kf.  There are 12
      % derivatives, associated with six qB and six qn DOFs.
      for idof = 1:6  % Start with the qB columns of Kf.

         icol = idref + idof;  % The column index of Kf.
         ic12 = 12*(idof-1);   % Indexes the derivative of qB.  Derivatives 1-6.

         % Start with the qB rows.
         rows = idref + [1:6].';
         dQhF(rows,icol) = dQhF(rows,icol) + (dQ(:,ic12+[1:6]).')*Fn;

         % Then the qn rows, if this is not the reference node.
         if (noddof ~= idref)
            rows = noddof + [1:6];
            dQhF(rows,icol) = dQhF(rows,icol) + (dQ(:,ic12+[7:12]).')*Fn;
         end

      end

      if (noddof ~= idref)

         for idof = 1:6  % Now the qn columns of Kf.

            icol = noddof + idof;  % The column index of Kf.
            ic12 = 12*(idof+6-1);  % Indexes the derivative of qB.  Derivatives 7-12.

            % The qB rows.
            rows = idref + [1:6].';
            dQhF(rows,icol) = dQhF(rows,icol) + (dQ(:,ic12+[1:6]).')*Fn;

            % Then the qn rows, if this is not the reference node.
            if (noddof ~= idref)
               rows = noddof + [1:6];
               dQhF(rows,icol) = dQhF(rows,icol) + (dQ(:,ic12+[7:12]).')*Fn;
            end

         end

      end

   end % for inod = 1:Nnod.

end % for ibod = 1:Nbod.

