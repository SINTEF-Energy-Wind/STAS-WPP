function df = getdf (f);

Nf = size(f,1);
df = zeros(Nf,1);
df(1) = f(2) - f(1);
df(Nf) = f(Nf) - f(Nf-1);
df(2:Nf-1) = 0.5*(f(3:Nf) - f(1:Nf-2));