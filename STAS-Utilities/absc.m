function z = absc(x)
ind = (real(x) < 0);
z = x;
z(ind) = -x(ind);
