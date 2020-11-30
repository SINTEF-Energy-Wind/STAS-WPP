function p = mvnorm (x,mu,sig)

r = chol(sig);
p = (2*pi)^(-size(x,1)/2) * exp (-sumsq (((x-mu).')*inv(r),2)/2) / prod (diag (r));