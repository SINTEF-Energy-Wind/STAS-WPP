function [dxdt,yout,A,B,C,D] = IBPLF (x,u,p,KpTab,KiTab)
%
% This function builds an IBP pitch controller which acts on measured
% blade root flap moments.  The moment measurements are low-pass
% filtered and MBC transformed to isolate the low-frequency asymmetric
% loading.  Filtering on the raw measurements used to derive the
% moments is not included as part of this function.
%
%   States:              y vector:             u vector:
%   Mpsim    1:2         bet      1:3          M      1:3
%   PsiM     3:4                               azi     4
%                                              sch     5
%                       
% Version:        Changes:
% --------        -------------
% 04.02.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 04.02.2019      Derivatives verified by complex step.
%
% Inputs:
% -------
% x               : 1,2: Mpsim (Nm)     Filtered blade root d-q moments.
%                   3,4: PsiM  (Nms)    Integral of Ki*Mpsim.
% u               : 1-3: M     (Nm)     Blade root moments.
%                   4:   azi   (rad)    Measured rotor azimuth.
%                   5:   W     (rad/s)  Measured rotor speed.
%                   6:   sch            Gain scheduling variable.
% p               : 1:   aLP   (rad/s)  LP filter frequency.

Nx = 4;
Nu = 5;
Ny = 3;

A  = sparse(Nx,Nx);
B  = sparse(Nx,Nu);
C  = sparse(Ny,Nx);
D  = sparse(Ny,Nu);

Mpsim = x(1:2);
PsiM  = x(3:4);
M     = u(1:3);
azi   = u(4);
sch   = u(5);
aLP   = p(1);

[Kp,dKp] = gains1 (sch,KpTab);
[Ki,dKi] = gains1 (sch,KiTab);

[TpsiB,TBpsi] = MBC (3,1,2,3,azi);
[dTpsiB,dTBpsi] = derivMBC (3,1,2,3,azi);

Mpsi = TBpsi(2:3,:)*M;  % cos and sin (d and q) components.

dxdt = [-aLP*Mpsim(1) + aLP*Mpsi(1); ...
        -aLP*Mpsim(2) + aLP*Mpsi(2); ...
        Ki*Mpsim(1);                 ...
        Ki*Mpsim(2)];

yout = TpsiB(:,2:3)*(Kp*Mpsim + PsiM);

A = [-aLP,    0, 0, 0; ...
        0, -aLP, 0, 0; ...
       Ki,    0, 0, 0; ...
        0,   Ki, 0, 0];

B = [aLP*TBpsi(2,1:3), aLP*dTBpsi(2,:)*M, 0; ...
     aLP*TBpsi(3,1:3), aLP*dTBpsi(3,:)*M, 0; ...
     0, 0, 0, 0, dKi*Mpsim(1);               ...
     0, 0, 0, 0, dKi*Mpsim(2)];


C = [TpsiB(:,2:3)*Kp, TpsiB(:,2:3)];

D = [zeros(3,3),                      ...
     dTpsiB(:,2:3)*(Kp*Mpsim + PsiM), ...
     TpsiB(:,2:3)*dKp*Mpsim];
  
