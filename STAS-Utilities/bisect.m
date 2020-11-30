function x = bisect (fun,x1,x2,del,maxit)

xlb = x1;
xub = x2;
ylb = fun(xlb);

x = 0.5*real(xlb + xub) + imag(xlb + xub);
y = fun(x);

iter = 0;
while ((abs(real(y)) > del) && (iter < maxit))

   iter = iter + 1;

   if (sign(real(y)) == sign(real(ylb)))
      % Move right.
      xlb = real(x);  % real: moving the boundary, so no longer retaining.
      ylb = y;        %       the complex part of input x1 or x2.
   else
      % Move left.
      xub = real(x);
   end

   x = 0.5*real(xlb + xub) + imag(xlb + xub);
   y = fun(x);

end
