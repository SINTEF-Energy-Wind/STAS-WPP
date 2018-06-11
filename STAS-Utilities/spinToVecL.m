function vec = spinToVecL (sp)

vec = zeros(3,1);
vec(1) =  sp(3,2);
vec(2) = -sp(3,1);
vec(3) =  sp(2,1);
