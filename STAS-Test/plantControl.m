function [dxdt,Phat,A,B,C,D] = plantControl (x,u,p,atables,cpct,linFlag)
%
% Build a wind plant control block associated with a given wind turbine.
%
% CAUTION: Complex step was tricky due to the numerical properties of the
% bisecting algorithm, in the solution for the target thrust.  The nominal
% calculation (for verification of A,B,C,D) works properly with complex
% step, but I have not checked all limiting cases, nor complex step with
% respect to parameters p.
%
%   States:              y vector:             u vector:
%   PsiP    1            Phat                  PhPCC    1
%   PsiT    2                                  PmPCC    2
%   uLP   3...24                               muV      3
%   n1    4...25                               Pa       4
%   n2    5...26                               FT       5
%                                              W        6
%                                              D        7
%                                              a1sum    8
%                                              a2sum    9
%
% Version:        Changes:
% --------        -------------
% 19.10.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 19.10.2020      A,B,C,D derivatives verified by complex step.
%
% Inputs:
% -------
% x               : 1: PsiP  (W)       Integral pathway on power control.
%                   2: PsiT  (W)       Integral pathway on load control.
% u               : 1: PhPCC (W)       Plant PCC power command.
%                   2: PmPCC (W)       Plant PCC measured power.
%                   3: muV   (m/s)     Cluster mean wind speed.
%                   4: Pa    (W)       The turbine's available power.
%                   5: FT    (N)       Observed rotor thrust (or other load).
%                   6: W     (rad/s)   Rotor speed.
%                   7: D     (-)       Damage accumulation rate.
%                   8: a1sum (-)       Sum of load alpha's from other turbines.
%                   9: a2sum (-)       Sum of dP alpha's from other turbines.
% p               : 1: KPP   (-)       KP on power control.
%                   2: KIP   (1/s)     KI on power control.
%                   3: KT    (m/s)     KT on load control.
%                   4: wlp   (rad/s)   Load branch LP filter corner frequency.
%                   5: PTLB  (W)       UB power deviation for load control. [not used]
%                   6: PTUB  (W)       LB power deviation for load control. [not used]
%                   7: weta  (-)       Width of the transition region in eta.
%                   8: Area  (m^2)     Rotor swept area.
%                   9: dens  (kg/m^3)  Air density.  Treated as a static parameter.
%                  10: alp   (rad/s)   Input LP filter corner frequency.
%                  11: an    (rad/s)   Notch filter frequency.
%                  12: zn1   (-)       Notch filter damping 1.  (zn1 << zn2)
%                  12: zn2   (-)       Notch filter damping 2.
% atables.a1,.a2  : alpha(D) functions.
% cpct            : Cp,Ct (V,W,b) table.
% linFlag         : Set to 1 to compute A,B,C,D.
%
% Outputs:
% --------
% dxdt            : Rate of change of states.
% Phat            : Turbine power command.
% A,B,C,D         : Linearized state matrices.

Nx = size(x,1);
Nu = size(u,1);
Ny = 1;

dxdt = zeros(Nx,1);
Phat = 0;

A = zeros(Nx,Nx);
B = zeros(Nx,Nu);
C = zeros(Ny,Nx);
D = zeros(Ny,Nu);

jlp = 3 + [0:3:24];
jn1 = 4 + [0:3:24];
jn2 = 5 + [0:3:24];

PsiP   = x(1);
PsiT   = x(2);
ulp    = x(jlp);
n1     = x(jn1);
n2     = x(jn2);
KPP    = p(1);
KIP    = p(2);
KT     = p(3);
wlp    = p(4);
PTLB   = p(5);
PTUB   = p(6);
weta   = p(7);
Area   = p(8);
dens   = p(9);
alp    = p(10);
an     = p(11);
zn1    = p(12);
zn2    = p(13);

a1table = atables.a1;
a2table = atables.a2;

% Filter all inputs.
an2  = an^2;
tz2a = 2*zn2*an;
tzza = 2*(zn1-zn2)*an;
dxdt(jlp) = -alp*ulp + alp*u;
dxdt(jn1) = n2;
dxdt(jn2) = -an2*n1 - tz2a*n2 + ulp;
ufilt = tzza*n2 + ulp;

PhPCC  = ufilt(1);
PmPCC  = ufilt(2);
muV    = ufilt(3);
Pa     = ufilt(4);
FT     = ufilt(5);
W      = ufilt(6);
Dam    = ufilt(7);
a1sum  = ufilt(8);
a2sum  = ufilt(9);

[alf1,da1dD] = gains1 (Dam,a1table);
lam1 = alf1/(alf1 + a1sum);

[alf2,da2dD] = gains1 (Dam,a2table);
lam2 = alf2/(alf2 + a2sum);

% Solve for the pitch angle that provides the desired power.
Pnom = lam1*PhPCC;  % The turbine's "ideal" share of the plant power.
CPnom = Pnom/(0.5*dens*Area*(muV^3));
RR = sqrt(Area/pi);
bmax = pi/6;




if (linFlag == 1)   % Temporary.





if (isreal ([CPnom, RR, muV, W]))  % No complex step.
   % First solve for the maximum Cp.
   [bmin,val,flg,outp] = fzero (@(b) getdcp (cpct,RR,muV,W,b), ...
                                [-pi/36, bmax]);
   [cpmax,ct,dcp,dct] = cpvwb (cpct,RR,muV,W,bmin);
   if (CPnom > cpmax) || (bmin >= bmax)
      % The desired point is cpmax or the min-pitch bound.
      bet = bmin;
      bminflg = true;
   else  % Solve for the pitch angle based on CPnom.
      [bet,val,flg,outp] = fzero (@(b) getcp (CPnom,cpct,RR,muV,W,b), ...
                                  [bmin, bmax]);
      bminflg = false;
   end
else
   % Custom routine suitable for complex step, but slower, thanks to 
   % Octave's interpreted-language overhead.
   bmin = bisect (@(b) getdcp (cpct,RR,muV,W,b),-pi/36,bmax,sqrt(eps),100);
   [cpmax,ct,dcp,dct] = cpvwb (cpct,RR,muV,W,bmin);
   if (real(CPnom) > real(cpmax)) || (real(bmin) >= real(bmax))
      bet = bmin;
      bminflg = true;
   else
      bet = bisect (@(b) getcp (CPnom,cpct,RR,muV,W,b),bmin,bmax,sqrt(eps),100);
      bminflg = false;
   end
end

[cp,ct,dcp,dct] = cpvwb (cpct,RR,muV,W,real(bet));  % Get final values and derivatives.
if (bminflg)
   % If the pitch angle is at the minimum boundary, then the gradient
   % with respect to CP is zero.
   dbetdCp = 0;
else
   dbetdCp = 1/dcp(3);
end
bet = real(bet) + i*(dbetdCp*(dcp(1)*imag(muV) + dcp(2)*imag(W)) ...
    +                (-dbetdCp)*imag(CPnom));  % Reconstruct the complex bet.
[cp,ct,dcp,dct] = cpvwb (cpct,RR,muV,W,bet);   % Complex step terms for later use.





else   % Temporary.





% TEMPORARY approximation to accelerate the calculation.
global betaglobal;
[cp,ct,dcp,dct] = cpvwb (cpct,RR,muV,W,betaglobal);
%deltacp = cp - CPnom;
%dbet = deltacp/dcp(3);  % dbeta = db/dcp dcp
%ct = ct + dct(3)*dbet;




end   % Temporary.






pAV2 = 0.5*dens*Area*(muV^2);
FhT = ct*pAV2;

epsP = lam2*(PhPCC - PmPCC);
epsT = FhT - FT;

dPP = KPP*epsP + PsiP;
dPT = PsiT;

Phat = lam1*PhPCC + dPP + dPT;

% starBlock has to be computed AFTER the output command.  It's OK,
% because the saturation doesn't take effect until the next timestep;
% it influences only dxdt, not y.
[etaP,detaP] = starBlock (Phat,0,Pa,epsP,weta,linFlag);

dxdt(1) = etaP*KIP*epsP;
dxdt(2) = -wlp*PsiT + wlp*KT*epsT;

if (linFlag == 1)

   dbetdmuV = dbetdCp*dcp(1);
   dbetdW   = dbetdCp*dcp(2);
   dbetdCPnom = -dbetdCp;

   aa1      =  alf1 + a1sum;
   dlam1dsa = -alf1/(aa1^2);
   dlam1da  =  (1/aa1) - alf1/(aa1^2);
   dlam1dD  =  dlam1da*da1dD;
   aa2      =  alf2 + a2sum;
   dlam2dsa = -alf2/(aa2^2);
   dlam2da  =  (1/aa2) - alf2/(aa2^2);
   dlam2dD  =  dlam2da*da2dD;

   dPnomdPhPCC = lam1;
   dPnomdsa    = PhPCC*dlam1dsa;
   dPnomdD     = PhPCC*dlam1dD;

   dFhTdW     = pAV2*(dct(2) + dct(3)*dbetdW);
   dFhTdmuV   = ct*dens*Area*muV                                   ...
              + pAV2*(dct(1)                                       ...
              +       dct(3)*(dbetdCPnom*(Pnom/pAV2)*(-3*muV^(-2)) ...
              +               dbetdmuV));

   dFhTdPnom  = pAV2*dct(3)*dbetdCPnom/(pAV2*muV);
   dFhTdPhPCC = dFhTdPnom*dPnomdPhPCC;
   dFhTdsa    = dFhTdPnom*dPnomdsa;
   dFhTdD     = dFhTdPnom*dPnomdD;

   depsTdFT    = -1;
   depsTdFhT   =  1;
   depsTdW     =  depsTdFhT*dFhTdW;
   depsTdmuV   =  depsTdFhT*dFhTdmuV;
   depsTdD     =  depsTdFhT*dFhTdD;
   depsTdsa    =  depsTdFhT*dFhTdsa;
   depsTdPhPCC =  depsTdFhT*dFhTdPhPCC;

   depsPdlam   =  PhPCC - PmPCC;
   depsPdD     =  depsPdlam*dlam2dD;
   depsPdsa    =  depsPdlam*dlam2dsa;
   depsPdPhPCC =  lam2;
   depsPdPmPCC = -lam2;

   dPhdepsP  =  KPP;
   dPhdlam1  =  PhPCC;
   dPhdD     =  dPhdlam1*dlam1dD +  dPhdepsP*depsPdD;
   dPhdsa1   =  dPhdlam1*dlam1dsa;
   dPhdsa2   =  dPhdepsP*depsPdsa;
   dPhdPhPCC =  lam1 + dPhdepsP*depsPdPhPCC;
   dPhdPmPCC =  dPhdepsP*depsPdPmPCC;

   detaPdPa    =  detaP(3);
   detaPdPh    =  detaP(1);
   detaPdepsP  =  detaP(4);

   Du = [dPhdPhPCC, dPhdPmPCC, 0, 0, 0, 0, dPhdD, dPhdsa1, dPhdsa2];

   Bu1 = (etaP*KIP + KIP*epsP*detaPdepsP) ...
       * [depsPdPhPCC, depsPdPmPCC, 0, 0, 0, 0, depsPdD, 0, depsPdsa] ...
       + KIP*epsP*(detaPdPh*Du + [0, 0, 0, detaPdPa, 0, 0, 0, 0, 0]);
   Bu2 = wlp*KT ...
       * [depsTdPhPCC, 0, depsTdmuV, 0, depsTdFT, depsTdW, depsTdD, depsTdsa, 0];

   A(1,1) = KIP*epsP*detaPdPh;
   A(1,2) = KIP*epsP*detaPdPh;
   A(2,2) = -wlp;
   
   A(jlp,jlp) = -alp*speye(Nu);     % dulp/dt = -alp*ulp + alp*u.
   A(jn1,jn2) = speye(Nu);          % dn1/dt = n2.
   A(jn2,jlp) = speye(Nu);          % dn2/dt = -an2*n1 - tz2a*n2 + ulp.
   A(jn2,jn1) = -an2*speye(Nu);
   A(jn2,jn2) = -tz2a*speye(Nu);
   A(1,jlp) = Bu1;                  % ufilt = tzza*n2 + ulp.
   A(1,jn2) = tzza*Bu1;
   A(2,jlp) = Bu2;
   A(2,jn2) = tzza*Bu2;

   B(jlp,:) = alp*speye(Nu);

   C(1,1:2) = [1, 1];
   C(1,jlp) = Du;
   C(1,jn2) = tzza*Du;

end % linFlag

end % plantControl




function f = getcp (CPnom,cpct,R,muV,W,bet)
   [cp,~,~,~] = cpvwb (cpct,R,muV,W,bet);
   f = cp - CPnom;
end


function dcp = getdcp (cpct,R,muV,W,bet)
   [~,~,dcps,~] = cpvwb (cpct,R,muV,W,bet);
   dcp = dcps(3);
end % getcp 