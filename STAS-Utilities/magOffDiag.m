function mofd = magOffDiag (T)
% Get the magnitudes of each skew-symmetric off-diagonal pair.

mofd = zeros(3,1);
mofd(1) = absc(T(3,2)) + absc(T(2,3));
mofd(2) = absc(T(1,3)) + absc(T(3,1));
mofd(3) = absc(T(1,2)) + absc(T(2,1));
