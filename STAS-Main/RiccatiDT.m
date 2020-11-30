function GG = RiccatiDT (meth,AA,BB,CC,QQ,RR,G0,tol)
%
% Integrates the matrix Riccati equation in time, or applies an 
% iterative solution, to find the steady-state solution.  Method
% 2 appears to be more reliable for discrete-time.
%
% Version:        Changes:
% --------        -------------
% 06.10.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 06.10.2020      
%
% Inputs:
% -------
% meth            : 1 = time integration.
%                   2 = iterative solution.
% AA,BB,CC        : State matrices.
% QQ,RR           : Process and measurement noise covariance matrices.
% G0              : Initial gain matrix.  (I - G0 C)A must be stable.
%
% Outputs:
% --------
% GG              : The gain matrix.  Nx-by-Ny.

itmax = 50;
iter = 0;
conv = 0;

Nx = size(AA,1);
Ny = size(CC,1);

GG = G0;
SS = speye(Nx);

if (meth == 1)

   while (conv == 0) && (iter <= itmax)
      iter = iter + 1;
      PH = (speye(Nx) - GG*CC)*AA;
[slap,shp,ifrq] = eigVal_silent(PH);
abs(slap).'
      SI = (speye(Nx) - GG*CC)*BB;
      SS1 = SS;
      SS = PH*SS1*(PH.') + SI*QQ*(SI.') + GG*RR*(GG.');
      ABC = (AA*SS*(AA.') + BB*QQ*(BB.'))*(CC.');
      GG = ABC*inv(CC*ABC + RR);
      del = max(max(abs(SS - SS1)));
      if (del < tol)
         conv = 1;
      end
printf('del = %+5.6e\n',del);
fflush(stdout);
   end

printf('conv = %4d\n',conv);

elseif (meth == 2)

   while (conv == 0) && (iter <= itmax)
      iter = iter + 1;
      SS1 = SS;
      PH = (speye(Nx) - GG*CC)*AA;
      SI = (speye(Nx) - GG*CC)*BB;
      RHS = SI*QQ*(SI.') + GG*RR*(GG.');
      SS = dlyap (PH,RHS);
      ABC = (AA*SS*(AA.') + BB*QQ*(BB.'))*(CC.');
      GG = ABC*inv(CC*ABC + RR);
      del = max(max(abs(SS - SS1)));
      if (del < tol)
         conv = 1;
      end
printf('del = %+5.6e\n',del);
fflush(stdout);
   end

[slap,shp,ifrq] = eigVal_silent(PH);
abs(slap).'
printf('conv = %4d\n',conv);

end


