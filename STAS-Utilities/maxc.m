function z = maxc (a,b)

ind = (real(a) < real(b));
z = a;
z(ind) = b(ind);