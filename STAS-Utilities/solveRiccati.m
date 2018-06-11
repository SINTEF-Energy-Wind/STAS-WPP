function PP = solveRiccati (BCflag,AA,BC,QQ,RR,TT)


% CAUTION, NOT ADAPTED FOR COMPLEX STEP DERIVATIVES.



%
% Integrates the matrix Riccati equation in time to find the steady-
% state solution.  A "lazy" approach is used, integrating the time-
% dependenent equations to the steady state.  This is functional but
% not the most efficient way to do it.
%
% Version:        Changes:
% --------        -------------
% 27.02.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 27.02.2017
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

solver = 2;  % 2 is much faster.  1 has been observed to be more
             % reliable when the weightings lead to poor numerical
             % conditioning.

if (solver == 1)

   global AA_Ric BC_Ric QQ_Ric RR_Ric

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

   ts = TT*[0 0.9 0.95 1]';
   Nt = size(ts,1);

   pvec0 = zeros(Nx^2,1);
   [pvec,istate,msg] = lsode (@RiccatiQS,pvec0,ts);

   'Infinite horizon Riccati equation'
%   ts'
%   pvec'

   PP = zeros(Nx,Nx);
   for icol = 1:Nx
      ind = Nx*(icol-1);
      PP(:,icol) = pvec(Nt,ind+[1:Nx]);
   end

%   clear AA_Ric BC_Ric QQ_Ric RR_Ric

elseif (solver == 2)

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

   [PP,LL,GG] = care (AA_Ric,BC_Ric,QQ_Ric,RR_Ric);

end

PP

'Riccati, should be zero'
Ric = (AA_Ric.')*PP + PP*AA_Ric + QQ_Ric ...
    - PP*BC_Ric*inv(RR_Ric)*(BC_Ric.')*PP;
'Normalized accuracy'
Ric./PP

clear AA_Ric BC_Ric QQ_Ric RR_Ric