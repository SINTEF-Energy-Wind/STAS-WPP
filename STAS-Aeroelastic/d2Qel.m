function d2Qe = d2Qel (d2Qu1,d2Qu2)
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
% d2Qu1,2         : A 6-by-12*12*12 matrix containing d2Qu/dq2 for each
%                   node.
%
% Outputs:
% --------
% d2Qe            : A 12-by-18*18*18 matrix containing d2Qe/dq2 for each 
%                   element.

d2Qe = zeros(12,18*18*18);

% Ref node - ref node.
for ii = 1:6
   ic324 = 324*(ii-1);      % Ref node DOFs.
   ic144 = 144*(ii-1);
   for jj = 1:6
      jc18 = 18*(jj-1);     % Ref node DOFs.
      jc12 = 12*(jj-1);
      d2Qe(:,ic324+jc18+[1:18]) = Qel(d2Qu1(:,ic144+jc12+[1:12]), ...
                                      d2Qu2(:,ic144+jc12+[1:12]));
   end
end

% Ref node - q1.
for ii = 1:6
   ic324 = 324*(ii-1);      % Ref node DOFs.
   ic144 = 144*(ii-1);
   for jj = 1:6
      jc18 = 18*(jj+6-1);   % q1 DOFs.
      jc12 = 12*(jj+6-1);
      d2Qe(:,ic324+jc18+[1:18]) = Qel(d2Qu1(:,ic144+jc12+[1:12]), ...
                                      zeros(6,12));
   end
end

% Ref node - q2.
for ii = 1:6
   ic324 = 324*(ii-1);      % Ref node DOFs.
   ic144 = 144*(ii-1);
   for jj = 1:6
      jc18 = 18*(jj+12-1);  % q2 DOFs.
      jc12 = 12*(jj+6-1);
      d2Qe(:,ic324+jc18+[1:18]) = Qel(zeros(6,12), ...
                                      d2Qu2(:,ic144+jc12+[1:12]));
   end
end

% q1 - Ref node.
for ii = 1:6
   ic324 = 324*(ii+6-1);    % q1 DOFs.
   ic144 = 144*(ii+6-1);
   for jj = 1:6
      jc18 = 18*(jj-1);     % Ref node DOFs.
      jc12 = 12*(jj-1);
      d2Qe(:,ic324+jc18+[1:18]) = Qel(d2Qu1(:,ic144+jc12+[1:12]), ...
                                      zeros(6,12));
   end
end

% q1 - q1.
for ii = 1:6
   ic324 = 324*(ii+6-1);    % q1 DOFs.
   ic144 = 144*(ii+6-1);
   for jj = 1:6
      jc18 = 18*(jj+6-1);   % q1 DOFs.
      jc12 = 12*(jj+6-1);
      d2Qe(:,ic324+jc18+[1:18]) = Qel(d2Qu1(:,ic144+jc12+[1:12]), ...
                                      zeros(6,12));
   end
end

% q1 - q2.  Zero.

% q2 - Ref node.
for ii = 1:6
   ic324 = 324*(ii+12-1);   % q2 DOFs.
   ic144 = 144*(ii+6-1);
   for jj = 1:6
      jc18 = 18*(jj-1);     % Ref node DOFs.
      jc12 = 12*(jj-1);
      d2Qe(:,ic324+jc18+[1:18]) = Qel(zeros(6,12), ...
                                      d2Qu2(:,ic144+jc12+[1:12]));
   end
end

% q2 - q1.  Zero.

% q2 - q2.
for ii = 1:6
   ic324 = 324*(ii+12-1);   % q2 DOFs.
   ic144 = 144*(ii+6-1);
   for jj = 1:6
      jc18 = 18*(jj+12-1);  % q2 DOFs.
      jc12 = 12*(jj+6-1);
      d2Qe(:,ic324+jc18+[1:18]) = Qel(zeros(6,12), ...
                                      d2Qu2(:,ic144+jc12+[1:12]));
   end
end

d2Qe = sparse(d2Qe);

