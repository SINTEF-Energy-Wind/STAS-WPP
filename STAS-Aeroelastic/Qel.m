function Qe = Qel (Qu1,Qu2)
%
% Make an element matrix from two nodal matrices.  Each nodal matrix
% of the type [v;w] = Qu (dq/dt) is 6-by-12, the 12 representing 6 body
% reference node DOFs and 6 nodal DOFs.  The element matrix is of the
% form
%              |O  |
%  |v1|        |Phi|
%  |w1| = [Qe] |d1 |
%  |v2|        |th1|
%  |w2|        |d2 |
%              |th2|
%
% and has dimensions 12-by-18.
%

Qe = zeros(12,18);

Qe(1:6,1:6)     = Qu1(:,1:6);
Qe(7:12,1:6)    = Qu2(:,1:6);
Qe(1:6,7:12)    = Qu1(:,7:12);
Qe(7:12,13:18)  = Qu2(:,7:12);
