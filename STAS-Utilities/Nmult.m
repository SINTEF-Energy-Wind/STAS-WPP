function b = Nmult (A,x)
%
% Multiply b = A*x, where A is blockwise diagonal consisting of uniform
% block sizes, such that the operation A*x involves m multiples of an
% N-by-N submatrix and a length N segment of x.
%
% (Trials indicate that this is slower than the nominal looping through
% m groups of matrix-vector multiplications.)
%
% Version:        Changes:
% --------        -------------
% 24.02.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 24.02.2018      
%
% Inputs:
% -------
% A               : Packed as N rows by m*N columns.
% x               : A vector of length m*N.
%
% Outputs:
% --------
% b               : Result of the multiplication.

N = size(A,1);
m = size(x,1)/N;

Adx = spalloc (size(A,1),size(A,2),nnz(A));
for icomp = 1:N
   Adx(icomp,:) = A(icomp,:).*(x.');
end

b = zeros(size(x,1),1);
for icomp = 1:N
   ics = icomp + [0:N:N*(m-1)].';
   vec = reshape(Adx(:,ics),size(x,1),1);
   b = b + vec;
end