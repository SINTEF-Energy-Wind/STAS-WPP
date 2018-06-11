function [Mp,rr,cr] = partitionMatrix (M,rDOFs,cDOFs)
%
% This function reorders the rows and columns of a matrix
% such that a specified list of DOFs is placed at the end.
%
% Version:        Changes:
% --------        -------------
% 11.12.2014      Original code.
%
% Version:        Verification:
% --------        -------------
% 11.12.2014      Verified on simple test cases.
%
% Inputs:
% -------
% M              : The matrix.
% rDOFs,cDOFs    : A vector of row and column DOFs to place
%                  at the end.  (Square matrix, identical
%                  reordering, set rDOFs=cDOFs.)  They must
%                  be in increasing order.
%
% Outputs:
% --------
% Mp             : The partitioned matrix.  Sparse.
% rr,cr          : Row and column ID lists of the DOFs which
%                  are retained.

% Here we can take advantage of the convention that
% M(vec1,vec2) returns a matrix consisting of the
% outer product of the rows listed in vec1 and the
% columns of vec2.
[Nr,Nc] = size(M);
Mp = sparse(Nr,Nc);

Nrdof = length(rDOFs);
Ncdof = length(cDOFs);

rr = zeros(Nr-Nrdof,1);
cr = zeros(Nc-Ncdof,1);
Nrr = Nr - Nrdof;
Ncr = Nc - Ncdof;

j = 1;
k = 1;
for idof = 1:Nr
   if (j <= Nrdof)
      if (idof == rDOFs(j))
         j++;
      else
         % Add to the list of non-reordered rows.
         rr(k) = idof;
         k++;
      end
   else
      % Add to the list of non-reordered rows.
      rr(k) = idof;
      k++;
   end
end

j = 1;
k = 1;
for idof = 1:Nc
   if (j <= Ncdof)
      if (idof == cDOFs(j))
         j++;
      else
         % Add to the list of non-reordered rows.
         cr(k) = idof;
         k++;
      end
   else
      % Add to the list of non-reordered rows.
      cr(k) = idof;
      k++;
   end
end

Mp(1:Nrr,1:Ncr) = M(rr,cr);
Mp(1:Nrr,Ncr+1:Nc) = M(rr,cDOFs);
Mp(Nrr+1:Nr,1:Ncr) = M(rDOFs,cr);
Mp(Nrr+1:Nr,Ncr+1:Nc) = M(rDOFs,cDOFs);



