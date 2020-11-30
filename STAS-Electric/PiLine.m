function [dxdt,A,By] = PiLine (Linflag,x,yin,params)
%
% The state equations for a Pi segment of an AC transmission line.
%
%   states:           y vector:
%   iLd,q    1,2      i(k+1)d,q        in   (Next cable segment)
%   vkd,q    3,4      v(k-1)d,q        in   (Previous cable segment)
%                     we               in   (Reference synchronous gen.)
%
% Version:        Changes:
% --------        -------------
% 29.08.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 29.08.2018      Derivatives verified by complex step.
%
% Inputs:
% -------
% Linflag         : Set to 1 for linearized equations.
% x               : 1,2: iL    (A)     Line currents
%                   3,4: vk    (V)     Node k (RHS) voltage
% yin             : 1,2: i_k+1 (A)     Next line segment current
%                   3,4: v_k-1 (V)     Node k-1 (LHS) voltage
%                   5:   we    (rad/s) Grid electrical frequency
% params          : 1:   R     (Ohms)  Line resistance
%                   2:   L     (H)     Line inductance
%                   3:   C     (F)     Node k capacitance
%
% Outputs:
% --------
%

Nx = 4;
Ny = 5;

A = sparse(Nx,Nx);
By = sparse(Nx,Ny);

R = params(1);
L = params(2);
C = params(3);
we = yin(5);

RL = R/L;
L1 = 1/L;
C1 = 1/C;

spn = [0 -1;1 0];

dxdt = [-(RL*eye(2) + we*spn)*x(1:2) + (yin(3:4) - x(3:4))*L1; ...
        -we*spn*x(3:4) + (x(1:2) - yin(1:2))*C1]; 

if (Linflag == 1)

   A(1:2,1:2)  = -RL*eye(2) - we*spn;
   A(1:2,3:4)  = -L1*eye(2);
   By(1:2,3:4) =  L1*eye(2);
   By(1:2,5)   = -spn*x(1:2);

   A(3:4,3:4)  = -we*spn;
   A(3:4,1:2)  =  C1*eye(2);
   By(3:4,1:2) = -C1*eye(2);
   By(3:4,5)   = -spn*x(3:4);

end 
