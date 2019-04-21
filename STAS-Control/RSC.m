function [dxdt,yout,A,B,C,D] = RSC (x,u,p,WVTab,WPTab,bminTab,KTables,mflag)
%
% Implement a rotor speed control function across the operating
% envelope.
%
%   States:              y vector:             u vector:
%   Wm      1            bhat                  W
%   betm    2            Pehat                 bet
%   PsiWb   3                                  V
%   Pf      4                                  Pc
%   Pehat   5
%   ntw    6,7
%   ntb    8,9
%   blp    10
%
% Version:        Changes:
% --------        -------------
% 28.01.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 28.01.2019      Derivatives verified by finite difference.
%
% Inputs:
% -------
% 
%
% Outputs:
% --------
% 

Nx = size(x,1);
Ny = 2;
Nu = size(u,1);

A = zeros(Nx,Nx);
B = zeros(Nx,Nu);
C = zeros(Ny,Nx);
D = zeros(Ny,Nu);

Wm    = x(1);
betm  = x(2);
PsiWb = x(3);
Pf    = x(4);
Pehat = x(5);
ntw   = x(6:7);
ntb   = x(8:9);
blp   = x(10);
W     = u(1);
bet   = u(2);
V     = u(3);
Pc    = u(4);
aW    = p(1);
ab    = p(2);
a1    = p(3);
a2    = p(4);
Pr    = p(5);
fcd   = p(6);
Wr    = p(7);
anw   = p(8);
z1nw  = p(9);
z2nw  = p(10);
anb   = p(11);
z1nb  = p(12);
z2nb  = p(13);
ablp  = p(14);

Wcd   = fcd*Wr;

aw2   = anw^2;
tz2aw = 2*z2nw*anw;
tzzaw = 2*(z1nw-z2nw)*anw;
ab2   = anb^2;
tz2ab = 2*z2nb*anb;
tzzab = 2*(z1nb-z2nb)*anb;

KpbTab = KTables(1);
KibTab = KTables(2);
KdbTab = KTables(3);

dxdt = zeros(Nx,1);

dxdt(1) = -aW*Wm + aW*W;
dxdt(2) = -ab*betm + ab*bet;

A(1,:) = [-aW, 0, 0, 0, 0, 0, 0, 0, 0, 0];
A(2,:) = [0, -ab, 0, 0, 0, 0, 0, 0, 0, 0];

B(1,:) = [aW, 0, 0, 0];
B(2,:) = [0, ab, 0, 0];

dxdt(6) = ntw(2);
dxdt(7) = -aw2*ntw(1) - tz2aw*ntw(2) + Wm;
dxdt(8) = ntb(2);
dxdt(9) = -ab2*ntb(1) - tz2ab*ntb(2) + betm;

Wfilt = tzzaw*ntw(2) + Wm;
dWfdt = tzzaw*dxdt(7) + dxdt(1);
bfilt = tzzab*ntb(2) + betm;
dbfdt = tzzab*dxdt(9) + dxdt(2);

dWfdx = [1, 0, 0, 0, 0, 0, tzzaw, 0, 0, 0];
ddWfdx = [-aW+tzzaw, 0, 0, 0, 0, -tzzaw*aw2, -tzzaw*tz2aw, 0, 0, 0];
ddWfdu = [aW, 0, 0, 0];

dbfdx = [0, 1, 0, 0, 0, 0, 0, 0, tzzab, 0];
ddbfdx = [0, -ab+tzzab, 0, 0, 0, 0, 0, -tzzab*ab2, -tzzab*tz2ab, 0];
ddbfdu = [0, ab, 0, 0];

A(6,:) = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0];
A(7,:) = [1, 0, 0, 0, 0, -aw2, -tz2aw, 0, 0, 0];
A(8,:) = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0];
A(9,:) = [0, 1, 0, 0, 0, 0, 0, -ab2, -tz2ab, 0];

[bmin,dbdV] = gains1 (V,bminTab);   % Maybe VLP?
bmax = pi/2;
bcen = 0.5*(bmax + bmin);

[Kpb,dKpb]  = gains1 (bfilt,KpbTab);
[Kib,dKib]  = gains1 (bfilt,KibTab);
[Kdb,dKdb]  = gains1 (bfilt,KdbTab);

% Solve for the transition speed as a function of the power, based
% on the curve of power as a function of speed.
[Wt,dWtdP] = solveNewt (@(W) gains1 (W,WPTab),Pc,0.8,eps,50,ones(50,1),50);

% Find the target speed for curtailed operation as a function of the
% observed windspeed.
if (real(Wt) >= real(Wcd))
   What = Wr;
   dWhdV = 0;
   dWhdP = 0;
else
   [What0,dWh0dV] = gains1 (V,WVTab);
   Wmax = 2*Wr;      % Arbitrary, should be > Wr by some reasonable amount.
   Wmin = Wt/fcd;    % A bit over the trans. speed, and by the definition Wcd = fcd*Wr
                     % this will be continuous with the other part of the "if".
   Wcen = 0.5*(Wmax + Wmin);
   WmWc = Wmax - Wcen;
   Wnorm = (What0 - Wcen)/WmWc;
   [Wns,dWns,d2Wns] = saturate (Wnorm,[0.99;1.01]);
   What = WmWc*Wns + Wcen;
   dWndWc = -(1/WmWc) + (What0 - Wcen)/(WmWc^2);
   dWhdV = dWns*dWh0dV;
   dWhdP = 0.5*WmWc*dWns*dWndWc*dWtdP/fcd ...
         + 0.5*(1 - Wns)*dWtdP/fcd;
end

% [From a perspective of Newton's method solution for operating points, 
% optimization, optimal control, etc. it would be much preferred to have
% a smooth transition rather than "if" statements...]

if ((((real(Wfilt) > real(Wcd)) || (real(Wfilt) > real(Wt))) && mflag == 0) || ...
    mflag == 1)
%printf('   RSC: Prescribed power and pitch control\n');
%fflush(stdout);
   % Prescribed power.
   Peh = Pc;

   % Pitch control.  Target speed is based on a windspeed estimate.
   epsW = Wfilt - What;
   bhat0 = Kpb*epsW + PsiWb + Kdb*dWfdt;

   dbhdx = Kpb*dWfdx + dKpb*epsW*dbfdx ...
         + [0, 0, 1, 0, 0, 0, 0, 0, 0, 0] ...
         + Kdb*ddWfdx + dKdb*dWfdt*dbfdx;
   dbhdu = [0, 0, -Kpb*dWhdV, -Kpb*dWhdP] ...
         + Kdb*ddWfdu;
   
   dxdt(4) = -a1*Pf + a1*Peh;
   A(4,:) = [0, 0, 0, -a1, 0, 0, 0, 0, 0, 0];
   B(4,:) = [0, 0, 0, a1];

   % Smooth saturation.
   bnorm = (bhat0 - bcen)/(bmax - bcen);
   [bns,dbns,d2bns] = saturate (bnorm,[0.99;1.01]);
   bhat = (bmax - bcen)*bns + bcen;

   % Anti-windup.
   bwn = (PsiWb - bcen)/(bmax - bcen);
   [bwns,dbwns,d2bwns] = saturate (bwn,[0.99;1.01]);   

   dxdt(3) = dbwns*Kib*epsW;
   A(3,:) = dbwns*Kib*dWfdx       ...
          + dbwns*dKib*epsW*dbfdx ...
          + d2bwns*Kib*epsW*[0, 0, 1, 0, 0, 0, 0, 0, 0, 0]/(bmax - bcen);
   B(3,:) = dbwns*[0, 0, -Kib*dWhdV, -Kib*dWhdP];

   dxdt(10) = -ablp*blp + ablp*bhat;
   A(10,:) = [0, 0, 0, 0, 0, 0, 0, 0, 0, -ablp] + ablp*dbns*dbhdx;
   B(10,:) = ablp*(dbns*dbhdu ...
           +      [0, 0, (dbns*(bhat0-bmax)*0.5*dbdV/(bmax - bcen))+(1-bns)*0.5*dbdV, 0]);
   

else
%printf('   RSC: Power from power-speed schedule\n');
%fflush(stdout);
   % Prescribed power from the operating schedule.
   [Peh,dPedW] = gains1 (Wfilt,WPTab);

   % Control pitch to the prescribed minimum angle.
   epsW = bmin - bfilt;
   bhat = Kpb*epsW + PsiWb - Kdb*dbfdt;

   dbhdx = -Kpb*dbfdx + dKpb*epsW*dbfdx ...
         + [0, 0, 1, 0, 0, 0, 0, 0, 0, 0] ...
         - Kdb*ddbfdx - dKdb*dbfdt*dbfdx;
   dbhdu = Kpb*[0, 0, dbdV, 0] ...
         - Kdb*ddbfdu;

   dxdt(3) = Kib*epsW;

   A(3,:) = -Kib*dbfdx + dKib*epsW*dbfdx;
   B(3,:) = [0, 0, Kib*dbdV, 0];

   dxdt(4) = -a1*Pf + a1*Peh;
   A(4,:) = [0, 0, 0, -a1, 0, 0, 0, 0, 0, 0] ...
          + a1*dPedW*dWfdx;

   dxdt(10) = -ablp*blp + ablp*bhat;
   A(10,:) = [0, 0, 0, 0, 0, 0, 0, 0, 0, -ablp] + ablp*dbhdx;
   B(10,:) = ablp*dbhdu;

end

dxdt(5) = -a2*Pehat + a2*Pf;

A(5,:) = [0, 0, 0, a2, -a2, 0, 0, 0, 0, 0];

C(1,:) = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1];
C(2,:) = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0];

yout = [blp;Pehat];

