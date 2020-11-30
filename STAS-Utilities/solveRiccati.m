function [PP,KK] = solveRiccati (solver,BCflag,AA,BC,QQ,RR,KK0,tf,dt)

% CAUTION, NOT ADAPTED FOR COMPLEX STEP DERIVATIVES.

%
% Integrates the matrix Riccati equation in time, or applies an 
% iterative solution, to find the steady-state solution.  The
% methods provided each have their strengths and weaknesses, so
% it is recommended to experiment and see which works best for
% a particular problem.
%
% Version:        Changes:
% --------        -------------
% 27.02.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 27.02.2017      Proven to drive the residual to zero, which by
%                 definition is the correct solution.
%
% Inputs:
% -------
% BCflag          : Set to 1 if this is a control problem and the
%                   matrix is B.  Set to 2 if this is an observer
%                   problem and the matrix is C.
% AA, BC          : The A matrix and either B or C matrix, with
%                   B for LQR and C for Kalman.
% QQ              : Weighting matrix on state vector, or "noise"
%                   covariance matrix forcing the state equation.
% RR              : Weighting matrix on the control, or "noise"
%                   covariance matrix on the measurements.
% TT              : Time interval to perform integration.
% 
%
% Outputs:
% --------
%

Nx = size(AA,1);

if (BCflag == 1)
   AA_Ric = AA;
   BC_Ric = BC;
   QQ_Ric = QQ;
   RR_Ric = RR;
elseif (BCflag == 2)
   AA_Ric = AA.';
   BC_Ric = BC.';
   QQ_Ric = QQ;   % Symmetric, regardless.
   RR_Ric = RR;   % Symmetric, regardless.
end

if (solver == 1)

   ts = [0:dt:tf].';
   Nt = size(ts,1);

   PP = zeros(Nx,Nx);
   for it = 1:Nt

if (mod(it,100) == 1)
printf ('%6d of %6d\n',it,Nt);
fflush(stdout);
end

      t = ts(it);
      [k1,jnk] = RiccatiQS (0,PP,AA_Ric,BC_Ric,QQ_Ric,RR_Ric);
      PP1 = PP + 0.5*k1*dt;
      [k2,jnk] = RiccatiQS (0,PP1,AA_Ric,BC_Ric,QQ_Ric,RR_Ric);
      PP = PP + k2*dt;

   end

elseif (solver == 2)

   % Caution, this tends to corrupt Octave's memory if the model
   % is large.  Problems were encountered with a few hundred DOFs.
   [PP,LL,GG] = care (AA_Ric,BC_Ric,QQ_Ric,RR_Ric);

elseif (solver == 3)

   flg = 0;
   iter = 0;
   itmax = 50;
   conv = sqrt(eps);
   KK = KK0;
   while ((flg == 0) && (iter <= itmax))
      iter = iter + 1;
      RHS = -(QQ_Ric + (KK.')*RR*KK);
      PP = lyap ((AA_Ric.')-(KK.')*(BC_Ric.'),-RHS);
      KK = RR\((BC_Ric.')*PP);
      Ric = (AA_Ric.')*PP + PP*AA_Ric - RHS;
      if (max(max(abs(Ric))) < conv)
         flg = 1;
      end
%printf('%5d  %+5.4e\n',iter,max(max(abs(Ric))));
   end

end

KK = RR\((BC_Ric.')*PP);

% Riccati residual, should be zero.
Ric = (AA_Ric.')*PP + PP*AA_Ric + QQ_Ric ...
    - PP*BC_Ric*inv(RR_Ric)*(BC_Ric.')*PP;

% Absolute accuracy, normalized accuracy.
[max(max(abs(Ric))) max(max(abs(Ric./PP)))]

