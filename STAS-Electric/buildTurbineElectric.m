function [dxdt,yout,A,B,C,D] = buildTurbineElectric (Linflag,x,yin,params)
%
% Build the unified state equations for a PMSG direct-drive generator with
% full power conversion, transformer, and converter controls.
%
%   states:           y vector:
% ------------------ External variables ---------------------
%                     wg         1     in
%                     we         2     in
%                     th_e       3     in
%                     vsd,q     4,5    in       yin
%                     ihgd,q    6,7    in
%                     Vhdc       8     in
%                     Qh         9     in
%                     Tg        10     out      yout
%                     isd,q    11,12   out
%
% ------------------- Local variables -----------------------
% ------------------------- PMSG --------------- 11 ---------
%   igd,q     1,2     vgd,q     1,2    in   (gen. current control)
%                     wg         3     in   (shaft)
%                     Tg         4     out
% ------ 2 ------- Converters and DC link ------ 15 ---------
%   Vdc        1      igd,q     1,2    in   (generator)
%                     vgd,q     3,4    in   (gen. current control)
%                     ipd,q     5,6    in   (transformer)
%                     vpd,q     7,8    in   (net. converter control)
%                     IIg        9     out
%                     IIn       10     out
% ------ 3 ------------- Transformer ----------- 25 ---------
%   ipd,q     1,2     vpd,q     1,2    in   (net. converter control)
%                     vsd,q     3,4    in   (grid)
%                     we         5     in   (grid)
%                     isd,q     6,7    out
% ------ 5 --------- Gen. current control ------ 32 ---------
%   imgd,q    1,2     wg         1     in   (shaft)
%   Psig      3,4     igd,q     2,3    in   (generator)
%   wemg       5      ihgd,q    4,5    in   (active power control)
%                     vgd,q     6,7    out
% ----- 10 ----------------- PLL --------------- 39 ---------
%   th_m       1      th_e       1     in   (grid)
%   vmsd,q    2,3     vsd,q     2,3    in   (grid)
%   Psie       4      wem        4     out
% ----- 14 -------- Reactive power control ----- 43 ---------
%   impd,q    1,2     ipd,q     1,2    in   (transformer)
%   imsd,q    3,4     isd,q     3,4    in   (transformer)
%   vmsd,q    5,6     vsd,q     5,6    in   (grid)
%   Vmdc       7      Vdc        7     in   (converter)
%   PsiDC      8      Vhdc       8     in   (input)
%   PsiQ       9      Qh         9     in   (voltage control)
%   PsiP     10,11    wem       10     in   (PLL)
%                     vpd,q    11,12   out
%
% Version:        Changes:
% --------        -------------
% 30.11.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 30.11.2018      Derivatives verified by complex step.
%
% Inputs:
% -------
% Linflag         : Set to 1 for linearized equations.
% x               : States at the current time, from list above.
% yin             : 1:    wg    (rad/s) Generator electrical speed.
%                   2:    we    (rad/s) Electrical speed.
%                   3:    th_e  (rad)   Electrical angle, integral of we.
%                   5-6:  vs    (V)     Network terminal voltage.
%                   7-8:  ihg   (A)     Commanded generator current.
%                   9:    Vhdc  (V)     Commanded DC link voltage.
%                   10:   Qh    (VA)    Commanded reactive power.
% params          : --------------------- PMSG ------------------------
%                   1:    np    (-)     Number of poles.
%                   2-5:  Ig    (A)     RMS phase currents for nonlinear inductance.
%                   6-9:  Lg    (H)     Nonlinear phase inductances.
%                   10:   Rg    (Ohms)  Phase resistance.
%                   11:   lamr  (Wb)    d-component flux linkage per phase.
%                   ------- 11 --------- DC link ----------------------
%                   1:    Cdc   (F)     DC link capacitance.
%                   2:    etac  (-)     Converter efficiency.
%                   ------- 13 ------- Transformer --------------------
%                   1:    a     (-)     Np/Ns turns ratio.
%                   2-5:  It    (A)     RMS phase currents for nonlinear inductance.
%                   6-9:  Lt    (H)     Lp + (a^2)Ls nonlinear phase inductances.
%                   10:   Rt    (Ohms)  Phase resistance.
%                   ------- 23 --- Gen. current control ---------------
%                   1-4:  Ig    (A)     RMS phase currents for nonlinear inductance.
%                   5-8:  Lg    (H)     Nonlinear phase inductances.
%                   9:    lamr  (Wb)    d-component flux linkage per phase.
%                   10:   KP    (1/s)   Proportional gain.
%                   11:   KI    (1/s^2) Integral gain.
%                   12:   KF    (-)     Feed-forward gain on flux linkage.
%                   13:   ag    (rad/s) Generator current filter frequency.
%                   14:   aw    (rad/s) Generator electric speed filter frequency.
%                   ------- 37 ------------ PLL -----------------------
%                   1:    av    (rad/s) Filter on voltage measurement.
%                   2:    KPe   (rad/s) Proportional gain.
%                   3:    KIe   (rad/s2) Integral gain.
%                   4:    weh   (rad/s) Reference frequency.
%                   ------- 41 -- Reactive power control --------------
%                   1:    ap    (rad/s) ip measurement filter.
%                   2:    ais   (rad/s) is measurement filter.
%                   3:    avs   (rad/s) vs measurement filter.
%                   4:    adc   (rad/s) Vdc measurement filter.
%                   5:    KpDC  (A/V)   Proportional gain on ipd control.
%                   6:    KiDC  (A/Vs)  Integral gain on ipd control.
%                   7:    KpQ   (1/V)   Proportional gain on ipq control.
%                   8:    KiQ   (1/Vs)  Integral gain on ipq control.
%                   9:    Kppd  (1/s)   Proportional gain on vpd control.
%                  10:    Kipd  (1/s^2) Integral gain on vpd control.
%                  11:    Kppq  (1/s)   Proportional gain on vpq control.
%                  12:    Kipq  (1/s^2) Integral gain on vpq control.
%                  13:    KFp   (-)     Feedforward gain on vms in vp control.
%                  14-17: It    (A)     RMS phase currents for nonlinear inductance.
%                  18-21: Lt    (H)     Lp + (a^2)Ls nonlinear phase inductances.
%
% Outputs:
% --------
% yout            : 1:    Tg    (Nm)    Generator air gap torque.
%                   2-3:  is    (A)     Network terminal currents.

Nxs  = [2 1 2 5 4 11].';
Nyis = [3 8 5 5 3 10].';
Nyos = [1 2 2 2 1  2].';
Nys  = Nyis + Nyos;
Nps  = [11 2 10 14 4 21].';
nms  = ['g';'d';'t';'p';'l';'r'];

Nx  = size(x,1);
Nyi = size(yin,1);                  % Size of yin.
Nyo = 3;                            % Size of yout.
Ny  = Nyi + Nyo;                    % Size of external [yin;yout].
Nyy = Ny + sum(Nyis) + sum(Nyos);   % Size of total: external plus local.

dxdt = zeros (Nx,1);
yout = zeros (Nyo,1);

% Define matrices for the internal un-linked y vector.
aa  = spalloc (Nx, Nx, 0.2*Nx*Nx);
bby = spalloc (Nx, Nyy,0.2*Nx*Nyy);
cc  = spalloc (Nyy,Nx, 0.2*Nx*Nyy);
ddy = spalloc (Nyy,Nyy,0.1*Nyy*Nyy);

% Indices.
ixg  = 0;
iye  = 0;
iyg  = Ny;
ipg  = 0;
for inm = 2:size(nms,1)
   eval(['ix'  nms(inm) ' = ix'  nms(inm-1) ' + Nxs(inm-1);']);
   eval(['iy'  nms(inm) ' = iy'  nms(inm-1) ' + Nyis(inm-1) + Nyos(inm-1);']);
   eval(['ip'  nms(inm) ' = ip'  nms(inm-1) ' + Nps(inm-1);']);
end
for inm = 1:size(nms,1)
   eval(['x'  nms(inm)       ' = x(ix'      nms(inm) ' + [1:Nxs(inm)]);']);
   eval([     nms(inm) 'params = params(ip' nms(inm) ' + [1:Nps(inm)]);']);
end

wg   = yin(1);
we   = yin(2);
th_e = yin(3);
vs   = yin(4:5);
ihg  = yin(6:7);
Vhdc = yin(8);
Qh   = yin(9);

ig   = xg;
ip   = xt;
Vdc  = xd;

% Build the individual components.  First the generator current control, as this
% depends only on the states x and external inputs yin.
ypin = [wg; ig; ihg];
[dxpdt,vg,Ap,Byp,Cp,Dyp] = converterControlGen (Linflag,xp,ypin,pparams);

% Then the PMSG.
ygin = [vg; wg];
[dxgdt,Tg,Ag,Byg,Cg]     = PMSG (Linflag,xg,ygin,gparams);

% Then the PLL.
ylin = [th_e; vs];
[dxldt,wem,Al,Byl,Cl]    = gridFrequency (Linflag,xl,ylin,lparams);

% Then the reactive power control.
is = tparams(1)*ip;
yrin = [ip; is; vs; Vdc; Vhdc; Qh; wem];
[dxrdt,vp,Ar,Byr,Cr,Dyr] = converterControlNet (Linflag,xr,yrin,rparams);

% Transformer.  The output 'is' has been computed above based on the transformer
% state variable.
ytin = [vp; vs; we];
[dxtdt,ist,At,Byt,Ct]  = Transformer (Linflag,xt,ytin,tparams);

% Finally the DC link.
ydin = [ig; vg; ip; vp];
[dxddt,II,Ad,Byd,Cd,Dyd] = DCLink (Linflag,xd,ydin,dparams);

% Load into the dx/dt vector.
for inm = 1:size(nms,1)
   eval(['dxdt(ix' nms(inm) ' + [1:Nxs(inm)]) = dx' nms(inm) 'dt;']);
end

% External outputs: 
yout(1)   = Tg;
yout(2:3) = ist;

if (Linflag == 1)

   % Linearized matrices.  First assemble the un-linked matrices.
   for inm = 1:size(nms,1)
      eval(['aa(ix' nms(inm) '+[1:Nxs(inm)],ix' nms(inm) ...
            '+[1:Nxs(inm)]) = A' nms(inm) ';']);
      eval(['bby(ix' nms(inm) '+[1:Nxs(inm)],iy' nms(inm) ...
            '+[1:Nys(inm)]) = By' nms(inm) ';']);
      eval(['cc(iy' nms(inm) '+[1:Nys(inm)],ix' nms(inm) ...
            '+[1:Nxs(inm)]) = C' nms(inm) ';']);
   end

   ddy(iyd+[1:size(Dyd,1)],iyd+[1:size(Dyd,2)]) = Dyd;
   ddy(iyp+[1:size(Dyp,1)],iyp+[1:size(Dyp,2)]) = Dyp;
   ddy(iyr+[1:size(Dyr,1)],iyr+[1:size(Dyr,2)]) = Dyr;

   % Fill in the connections with the cc and ddy matrices.  
   ddy(iyg+[1:2],iyp+[6:7]) = speye(2);     % PMSG vg to gen. cont. vg.
   ddy(iyg+3,iye+1) = 1;                    % PMSG wg to external wg.
   cc(iyd+[1:2],ixg+[1:2]) = speye(2);      % DC link ig to PMSG ig.
   ddy(iyd+[3:4],iyp+[6:7]) = speye(2);     % DC link vg to gen. cont. vg.
   cc(iyd+[5:6],ixt+[1:2]) = speye(2);      % DC link ip to trans. ip.
   ddy(iyd+[7:8],iyr+[11:12]) = speye(2);   % DC link vp to Q cont. vp.
   ddy(iyt+[1:2],iyr+[11:12]) = speye(2);   % Trans. vp to Q cont. vp.
   ddy(iyt+[3:4],iye+[4:5]) = speye(2);     % Trans. vs to external vs.
   ddy(iyt+5,iye+2) = 1;                    % Trans. we to external we.
   ddy(iyp+1,iye+1) = 1;                    % gen. cont. wg to external wg.
   cc(iyp+[2:3],ixg+[1:2]) = speye(2);      % gen. cont. ig to PMSG ig.
   ddy(iyp+[4:5],iye+[6:7]) = speye(2);     % gen. cont. ihg to external ihg.
   ddy(iyl+1,iye+3) = 1;                    % PLL th_e to external th_e.
   ddy(iyl+[2:3],iye+[4:5]) = speye(2);     % PLL vs to external vs.
   cc(iyr+[1:2],ixt+[1:2]) = speye(2);      % Q cont. ip to trans. ip.
   ddy(iyr+[3:4],iyt+[6:7]) = speye(2);     % Q cont. is to trans. is.
   ddy(iyr+[5:6],iye+[4:5]) = speye(2);     % Q cont. vs to external vs.
   cc(iyr+7,ixd+1) = 1;                     % Q cont. Vdc to DC link Vdc.
   ddy(iyr+8,iye+8) = 1;                    % Q cont. Vhdc to external Vhdc.
   ddy(iyr+9,iye+9) = 1;                    % Q cont. Qh to external Qh.
   ddy(iyr+10,iyl+4) = 1;                   % Q cont. wem to PLL wem.

   % External outputs.
   ddy(iye+Nyi+1,iyg+4) = 1;                % Ext. Tg to PMSG Tg.
   ddy(iye+Nyi+[2:3],iyt+[6:7]) = speye(2); % Ext. is to trans. is.

   % Rearrange the matrices in terms of a 'u' vector consisting of the external
   % inputs, and the remaining 'y' vector.
   indu = [1:Nyi].';
   indy = [Nyi+1:Nyy].';
   bbu  = bby(:,indu);
   ddu  = ddy(indy,indu);
   bby  = bby(:,indy);
   ddy  = ddy(indy,indy);
   cc   = cc(indy,:);

   Dnorm = ones(Nyy-Nyi,1);
   [A,B,C,D] = modularToUnifiedStateSpace (aa,bbu,bby,cc,ddu,ddy,Dnorm);

else

   Nyr = Nyy - Nyi;
   A = sparse(Nx,Nx);
   B = sparse(Nx,Nyr);
   C = sparse(Nyr,Nx);
   D = sparse(Nyr,Nyr);

end

