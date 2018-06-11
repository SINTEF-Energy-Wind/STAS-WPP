function vec = spinToVec (sp)

vec = zeros(3,1);
vec(1) = 0.5*(sp(3,2) - sp(2,3));
vec(2) = 0.5*(sp(1,3) - sp(3,1));
vec(3) = 0.5*(sp(2,1) - sp(1,2));
