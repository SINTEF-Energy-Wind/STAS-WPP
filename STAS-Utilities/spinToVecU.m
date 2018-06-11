function vec = spinToVecU (sp)

vec = zeros(3,1);
vec(1) = -sp(2,3);
vec(2) =  sp(1,3);
vec(3) = -sp(1,2);
