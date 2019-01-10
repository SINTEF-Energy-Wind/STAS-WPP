function [dxdt,yout,A,By,C] = Transformer (LinFlag,x,yin,params)
%
% A simple impedance model of a transformer, with the possibility to
% define an effective nonlinear inductance to account for saturation,
% particularly for overpower operation.
%
%   states:           y vector:
%   ipd,q     1,2     vpd,q     1,2    in   (converter control)
%                     vsd,q     3,4    in   (grid)
%                     we         5     in   (grid)
%                     isd,q     6,7    out
%
% Version:        Changes:
% --------        -------------
% 05.09.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 05.09.2018      Derivatives verified with complex step.
%
% Inputs:
% -------
% Linflag         : Set to 1 for linearized equations.
% x               : 1,2: ip    (A)     Primary winding current
% yin             : 1,2: vp    (V)     Primary terminal voltage
%                   3,4: vs    (V)     Secondary terminal voltage
%                    5:  we    (rad/s) Electrical frequency
% params          :  1:  a     (-)     Np/Ns turns ratio
%                   2-5: It    (A)     RMS phase currents for nonlinear inductance
%                   6-9: Lt    (H)     Lp + (a^2)Ls nonlinear phase inductances
%                   10:  Rt    (Ohms)  Phase resistance
%
% Outputs:
% --------
% 

Nx = 2;
Ny = 7;

L  = sparse(Nx,Nx);
A  = sparse(Nx,Nx);
By = sparse(Nx,Ny);
C  = sparse(Ny,Nx);

ip = x;

vp = yin(1:2);
vs = yin(3:4);
we = yin(5);

a  = params(1);
It = params(2:5);
Lt = params(6:9);
Rt = params(10);

two = 2*pi/3;
four = 4*pi/3;
sq = sqrt(2/3);
sq2 = sq/sqrt(2);

psi = 0;  % dq is independent of this.
th = 0;

ct = cos(th);
st = sin(th);
c2 = cos(th-two);
s2 = sin(th-two);
c4 = cos(th-four);
s4 = sin(th-four);

lamr = params(2)*[ct;c2;c4];

Tat = sq*[ct  c2  c4; ...
         -st -s2 -s4];
Tta = sq*[ct -st; ...
          c2 -s2; ...
          c4 -s4];
dTtadth = sq*[-st -ct; ...
              -s2 -c2; ...
              -s4 -c4];

% Compute the RMS phase current and phase inductance.
magI = maxc(sqrt(ip(1)^2 + ip(2)^2),eps);
Irms = magI*sq2;
[Lph,dLph] = inductance (It,Lt,Irms);

Lambda  = Lph*eye(3);
dLambda = dLph*eye(3);
dLamN   = [dLambda*ip(1)*sq2/magI dLambda*ip(2)*sq2/magI];
TLT     = Tat*Lambda*Tta;
TLamDT  = Tat*Lambda*dTtadth;
WTLamDT = we*TLamDT;
RR      = Rt*eye(3);
TRT     = Tat*RR*Tta;

dxdt = TLT\(-(TRT + WTLamDT)*ip + vp - a*vs);
yout = a*ip;

L(:,:)    =  TLT;

A(:,:)    = -(TRT + WTLamDT);
vec       =  Tta*dxdt + we*dTtadth*ip;
A(:,1)    =  A(:,1) - Tat*dLamN(:,1:3)*vec;
A(:,2)    =  A(:,2) - Tat*dLamN(:,4:6)*vec;

By(:,1:2) =  eye(2);
By(:,3:4) = -a*eye(2);
By(:,5)   = -TLamDT*ip;

C(6:7,:)  =  a*eye(2);

A = L\A;
By = L\By;
