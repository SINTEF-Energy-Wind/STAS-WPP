function F = RCF (t,t0,t1,Fmax,Ndof)

   F = zeros(Ndof,1);
   if ((t > t0) && (t < t1))
      F(Ndof-6+2) = Fmax*(sin((t-t0)*pi/(t1-t0))^2);
   end
