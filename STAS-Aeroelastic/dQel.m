function dQe = dQel (dQu1,dQu2)
%
% Version:        Changes:
% --------        -------------
% 13.11.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 13.11.2017      
%
% Inputs:
% -------
% dQu1,2          : A 6-by-12*12 matrix containing dQu/dq for each node.
%
% Outputs:
% --------
% dQe             : A 12-by-18*18 matrix containing dQe/dq for each element.

dQe = zeros(12,18*18);

for ii = 1:6
   ic18 = 18*(ii-1);     % Ref node DOFs.
   ic12 = 12*(ii-1);
   dQe(:,ic18+[1:18]) = Qel(dQu1(:,ic12+[1:12]),dQu2(:,ic12+[1:12]));
end

for ii = 1:6
   ic18 = 18*(ii+6-1);   % Node 1 DOFs.
   ic12 = 12*(ii+6-1);
   dQe(:,ic18+[1:18]) = Qel(dQu1(:,ic12+[1:12]),zeros(6,12));
end

for ii = 1:6
   ic18 = 18*(ii+12-1);  % Node 2 DOFs.
   ic12 = 12*(ii+6-1);
   dQe(:,ic18+[1:18]) = Qel(zeros(6,12),dQu2(:,ic12+[1:12]));
end

dQe = sparse(dQe);

