% Run a Newton's method solution of the closed-loop turbine, and check the
% control tuning.  This attempts a single-level solution, that is, calling
% one function that returns dx/dt and A for the entire turbine.  It is also
% possible to pursue a hierarchical solution where each module is solved
% independently with inner Newton loops, and then the overall system with
% an outer Newton loop acting on the residuals in the interface variables.
%
% asymWind is set up to accommodate yaw angles and asymmetric wind fields,
% in addition to the nominal case of uniform aligned wind.

clear;

%============================================
% Abbreviated inputs for uniform wind.
ffllgg = 1;
if (ffllgg == 1)
   initFlagin = 1;
   lcin = '_P060_V100';
   outnmin = '_P060_';
   Vinp = 10.;
   RSCin = 1;    % 1: Pitch, 2: Variable-speed.
   Pcin = 6.e6;
   drin = 1;     % 1: full, 3: rigid struct, 5: elastic struct, 6: aero only.
elseif (ffllgg == 2)
   initFlagin = 1;
   lcin = '_V055';
   outnmin = '_';
   Vinp = 5.5;
   RSCin = 2;    % 1: Pitch, 2: Variable-speed.
   Pcin = 1.e7;
   drin = 1;     % 1: full, 3: rigid struct, 5: elastic struct, 6: aero only.
elseif (ffllgg == 3)
   initFlagin = 1;
   lcin = '_P070_V100';
   outnmin = '_P070_';
   Vinp = 10.0;
   RSCin = 1;    % 1: Pitch, 2: Variable-speed.
   Pcin = 7.e6;
   drin = 1;     % 1: full, 3: rigid struct, 5: elastic struct, 6: aero only.
end
%============================================

printf('Reading inputs\n');
fflush(stdout);

nm = 'DTU10MW';
initFlag = initFlagin;  % 0: Initial guess.  1: Load from file (using lcnm ID).
lcnm = lcin;
outnm = outnmin;

eval(['[s,a] = STASTurbine_'  nm ' ();']);
eval(['epar  = STASElectric_' nm ' ();']);
eval(['ppar  = STASPitch_'    nm ' ();']);
eval(['ypar  = STASYaw_'      nm ' ();']);
eval(['c     = STASControl_'  nm ' ();']);

[length,time,mass,current,voltage,        ...
 velocity,force,power,stress,ndens,nvisc, ...
 stiffness,damping,resistance,inductance, ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
[imdofs,Nmd] = getmdofRefs (s);

Nb  = a.Nb;
Neb = a.Neb;
Nae = Nb*Neb;

% For linear analysis, the azimuth angle is arbitrary, as its
% influence will be eliminated by the MBC transform.
azi0 = 0*pi/180;
cp0 = cos(azi0);
sp0 = sin(azi0);

% ===================================================================
% Input parameters defining the load case.
%
% Note that the tower mode shapes will be slightly more accurate if
% an initial yaw angle is input via yawr; then the mode shapes will
% be computed using this yaw position, rather than zero.
%
a.dens = 1.225                     / ndens;      % Air density.
a.visc = 1.789e-5                  / nvisc;      % Air viscosity.
Vmax   = Vinp                       / velocity;
Vmin   = Vinp                       / velocity;
Vmag   = 0.5*(Vmax + Vmin);
yawr   = 0*pi/180;                               % P: zero-q nodal disp.
betar  = 0*pi/180;                               % P: zero-q nodal disp.
yawg   = 0*pi/180;                               % q0 relative to bref.
betag  = 5*pi/180;                               % q0 relative to bref.
grav   = [0;0;-9.807]              / (length/(time^2));

c.RSCFlag = RSCin;  % Force the flag to correspond to the wind speed.
                % 1: pitch control.  2: variable-speed.

vs     = [33000;0]                 / voltage;    % Grid electrical voltage.
we     = 50*(2*pi)                 * time;       % Grid electrical frequency.
th_e   = 0;                                      % Ref. elec. angle.
Vhdc   = c.Vhdc;                                 % DC link voltage command.
Qh     = 0                         / power;      % Reactive power command.

Area   = pi*(c.Ro^2);
Pguess = minc(0.45*0.5*a.dens*Area*(Vmag^3),c.Pr);
Pc     = Pcin/power; % c.Pr;                                   % No power command.
%Pc     =  7.e6                     / power;      % Power command, use
%Pguess = Pc;                                     % these lines.
[Wguess,dW] = gains1 (Vmag,c.WVTab);

wg     = 0.5*c.np*Wguess;                        % Generator elec. speed.
ihg    = [1/current;-Pguess/(3000/voltage)];     % Gen. current commands.

betars = [betar;betar;betar];
betags = [betag;betag;betag];

% ===================================================================
% Solve for the initial electrical system states.
printf('Solving for initial electrical states\n');
fflush(stdout);

arat = epar(14);
ip   = [Pguess/(vs(1)*arat);1/current];
yein = [wg;we;th_e;vs;ihg;Vhdc;Qh];
xein = [ihg;Vhdc;vs*arat;ihg;ihg;we;th_e;vs;0; ...
        ip;ip*arat;vs;Vhdc;Vhdc;Qh;ip];

cnv = eps^0.6;
Ns = 50;
bta = ones(Ns,1);
litmax = 20;
[xe0,jnk] = solveNewt (@(x) BTEfun(x,yein,epar),0,xein,cnv,Ns,bta,litmax);

% ===================================================================
% Initialize the structral calculation or load initial values and
% mode shapes.
printf('Initializing structural states and mode shapes\n');
fflush(stdout);

[xs0,etas0,q0,dq0dt,d2q0dt2,P,shape0,freq0,mdamp0,ret,slv, ...
 Ndj,Nnod,Neta] = structInit (s,yawr,azi0,betars,Wguess);
Nret = size(ret,1);
Nslv = size(slv,1);
Nylin = 4*Ndj + 9*Nnod + 137*Nae + 5;
Nu = Ndj + 3*Nae;

% Put the guessed pitch angle into the eta vector.
etas0(Neta-[2 1 0]) = betags;

% ===================================================================
% BEM setup, parameters that stay fixed throughout the calculations,
% and initial estimates of the aerodynamic states.
printf('Initializing aerodynamic states\n');
fflush(stdout);

[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);
[Tas,ch,Lel,foilwt,aoaz,aoast,xas,yas,iq] = BEMsetup (s,a);

Td_n = [cp0 -sp0 0;sp0 cp0 0;0 0 1];
Try = Tn_y;
Tyy0 = TFromTheta (q0(idofs(3)+[4:6]));
Ty0g = TFromTheta (P(idofs(3)+[4:6]));
Trg = Ty0g*Tyy0*Try;

[Tar,Trg,TB0g,TsB,TBB0,dTar,dTsB,dTBB0,wga] = ...
                      BEMprepTransforms (s,a,q0,dq0dt,P,Tas);

[zr,Area,Dp,rp,Lp,xeg,xhg,xyg] = ...
         BEMprepProjections (s,a,q0,P,Try,Trg);
Dps = zeros(Nae,1);
Dps(1:Neb)         = Dp(1);
Dps(Neb+[1:Neb])   = Dp(2);
Dps(2*Neb+[1:Neb]) = Dp(3);

% ===================================================================
% Fill out the global velocity vector based on the inputs.
% Here the wind has an asymmetry.  These inputs are in body
% coordinates and are later transformed to MBC.
Vg = zeros(3*Nae,1);
thel = atan2(zr(2:2:2*Nae),zr(1:2:2*Nae-1));
cthel = cos(thel);
sthel = sin(thel);

% Horizontal gradient (partial wake state, for instance).
Vg(1:3:3*Nae-2) = 0.5*(Vmax + Vmin) + 0.5*(Vmax - Vmin)*rp.*cthel/rp(Neb);

% Vertical gradient (wind shear).
%Vg(1:3:3*Nae-2) = 0.5*(Vmax + Vmin) + 0.5*(Vmax - Vmin)*rp.*sthel/rp(Neb);

Vg(2:3:3*Nae-1) = 0;
Vg(3:3:3*Nae)   = 0;
%[thel cthel sthel Vg(1:3:3*Nae-2)]
% ===================================================================

% Initial guess for induced velocity.
Viguess = zeros(2*Nae,1);
Viguess(1:2:2*Nae-1) = -0.35*Vg(1:3:3*Nae-2);

xa0 = BEMinit (1,                                ...
               Viguess,Tar,Trg,Vg,wga,zr,ch,Lel, ...
               aoast,a.dens,Area,Dps,azi0,Wguess);

% Best estimate of the initial modal aero states.
bsh = bladeModeShape (s,ret,shape0);
Psi = aeroPsi (a,rp,bsh);
etaa0 = pinv(Psi)*xa0;
Naero = size(Psi,2);

% ===================================================================
% Actuators.
printf('Initializing actuator states\n');
fflush(stdout);

yp0 = zeros(9,1);
yp0([1 3 5]) = betags;
yp0([2 4 6]) = betags;
xp0 = zeros(6,1);
xp0([1 3 5]) = betags;
yy0 = zeros(3,1);
xy0 = zeros(2,1);

% ===================================================================
% Initialize the control states.
printf('Initializing control states\n');
fflush(stdout);

bavg = sum(betags)/3;

xc0 = [Wguess;Vmag;Wguess;bavg;bavg;Pguess;Pguess;Wguess;Wguess;bavg;bavg;bavg; ...
       0;0;Pguess;ihg(2);0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

% ===================================================================
% Initial conditions.
printf('Assembling initial conditions.\n');
fflush(stdout);

igen = [idofs(3)+[1:6] idofm(5)+[1:6] idofs(4)+[1:6] idofs(5)+[1:6]].';
ipit = [Ndof+[4:6]].';
iyaw = Ndof + 1;

u0 = [zeros(Ndj,1);zeros(Neta,1);Vg;we;th_e;vs;Pc;Qh];

% Initial states.
x0 = [etas0;etaa0;xp0;xy0;xe0;xc0];    % Guess
dx0 = zeros(size(x0,1),1);
dx0(1:Neta) = x0(Neta+[1:Neta]);

Nx = size(x0,1);
Nu = size(u0,1);

% ===================================================================
% MBC transforms.
printf('Building MBC transforms.\n');
fflush(stdout);

Nxa = size(Psi,2);
Nxp = 6;
[blxdof,bludof,blydof] =  MBCindices (Nxa,Nxp,Neta,Ndj,Neb,imdofs,idofs);
[TpsixB,TBxpsi]   = MBC      (Nx,blxdof(:,1),blxdof(:,2),blxdof(:,3),azi0);
[TpsiuB,TBupsi]   = MBC      (Nu,bludof(:,1),bludof(:,2),bludof(:,3),azi0);
bldof = {blxdof,bludof,blydof};

% ===================================================================
% Initialize the MBC variables.

iueta = Ndj + [1:Neta].';
niueta = [[1:Ndj] [Ndj+Neta+1:Nu]].';
upsi = TBupsi(niueta,:)*u0;

if (initFlag == 0)

printf('Initial conditions are guessed.\n');
fflush(stdout);

   xpsi = TBxpsi*x0;

   % The initial value of dxpsi should match xpsi, via the appropriate
   % transform, in the velocity DOFs.  To do this, convert to body
   % coordinates, assign dx, and transform back to dxpsi.
   iazi = Neta-4;
   [x,Txpx,jnk2] = buildTxpx (xpsi,zeros(Nx,1),blxdof(:,1),blxdof(:,2),blxdof(:,3),iazi);
   dx = zeros(Nx,1);
   dx(1:Neta) = x(Neta+[1:Neta]);
   dxpsi = Txpx\dx;

   % If desired, implement gravity as a constant external force.  This
   % may aid convergence of the structural DOFs.  Comment out if gravity
   % is to be considered as an integral part of the element matrices,
   % see structure.m.
   g0 = [0;0;0];
   [Lpsi,Rpsi,ypsi,Apsi,Bpsi,Cpsi,Dpsi] =                   ...
            MBCCLT (0,xpsi,dxpsi,upsi,s,a,                  ...
                    epar,ppar,ypar,c,g0,P,shape0,mdamp0,    ...
                    Tas,Try,ch,Lel,foilwt,aoaz,aoast,       ...
                    xas,yas,Psi,igen,ipit,iyaw);

   Rgrav = Lpsi(:,Neta+[1:3])*grav;  % Mass matrix columns times accel.
   grav = [0;0;0];

else

printf('Initial conditions are read from a file.\n');
fflush(stdout);

   eval(["load 'xpsi" lcnm ".txt';"]);
   eval(["load 'dxpsi" lcnm ".txt';"]);
   eval(["load 'Rgrav" lcnm ".txt';"]);
   eval(['xpsi  = xpsi' lcnm ';']);
   eval(['dxpsi = dxpsi' lcnm ';']);
   eval(['Rgrav = Rgrav' lcnm ';']);

end


drflg = drin;

if (drflg == 1)
   % Locked azimuth (free W).  Tower damping control functions, no VIG, IPB.
   dret = [[1:Neta-6] Neta-[5 3 2 1 0] Neta+[1:Neta-6] 2*Neta-[5 4 3 2 1 0] ...
           2*Neta+[1:Naero] 2*Neta+Naero+[1:8] 2*Neta+Naero+8+[1:25] ...
           2*Neta+Naero+8+25+[[1:12] [15:24]]].';

elseif (drflg == 2)
   % Locked azimuth (free W).  No tower damping, VIG, IPB.
   dret = [[1:Neta-6] Neta-[5 3 2 1 0] Neta+[1:Neta-6] 2*Neta-[5 4 3 2 1 0] ...
           2*Neta+[1:Naero] 2*Neta+Naero+[1:8] 2*Neta+Naero+8+[1:25] ...
           2*Neta+Naero+8+25+[[1:12] [15:18]]].';

elseif (drflg == 3)
   % Rigid structure.  It is good to start a new analysis with a rigid
   % solution, which allows the rotor speed, control, and electric parameters
   % to get close to their final values.
   dret = [Neta-[5 3 2 1 0] 2*Neta-[5 4 3 2 1 0] ...
           2*Neta+[1:Naero] 2*Neta+Naero+[1:8] 2*Neta+Naero+8+[1:25] ...
           2*Neta+Naero+8+25+[1 2 3 4 5 6 7 8 9 10 11 12 15 16 17 18]].';

elseif (drflg == 4)
   % No aero.
   dret = [[1:Neta-6] Neta-[5 3 2 1 0] Neta+[1:Neta-6] 2*Neta-[5 4 3 2 1 0] ...
           2*Neta+Naero+[1:8] 2*Neta+Naero+8+[1:25] ...
           2*Neta+Naero+8+25+[[1:12] [15:24]]].';

elseif (drflg == 5)
   % Only elastic structural.
   dret = [[1:Neta-6] Neta+[1:Neta-6]].';

elseif (drflg == 6)
   % Only aero.  If the optimization is "sticking" at a nonzero metric,
   % try running a pure-aero solution followed by an elastic structural
   % solution, then close out with the full solution.
   dret = [2*Neta+[1:Naero]].';

elseif (drflg == 7)

%   dret = [Neta-[5 3 2 1 0] 2*Neta-[5 4 3 2 1 0] ...
%           [209:288] [293:334] 2*Neta+Naero+[1:8] 2*Neta+Naero+8+[1:25] ...
%           2*Neta+Naero+8+25+[1 2 3 4 5 6 7 8 9 10 11 12 15 16 17 18]].';
%   dret = [[1:Neta-6] Neta-[5 3 2 1 0] Neta+[1:Neta-6] 2*Neta-[5 4 3 2 1 0] ...
%           [209:288] [293:334] 2*Neta+Naero+[1:8] 2*Neta+Naero+8+[1:25] ...
%           2*Neta+Naero+8+25+[[1:12] [15:24]]].';
%   dret = [[209:288] [293:334]].';
%   dret = [291 292].';


end

Ndr = size(dret,1);

% ===================================================================
% Get the initial N-R residual.
printf('Newton-Raphson solution.\n');
printf('   Computing the initial residual.\n');
fflush(stdout);

[Lpsi,Rpsi,ypsi,Apsi,Bpsi,Cpsi,Dpsi] =                   ...
         MBCCLT (0,xpsi,dxpsi,upsi,s,a,                  ...
                 epar,ppar,ypar,c,grav,P,shape0,mdamp0,  ...
                 Tas,Try,ch,Lel,foilwt,aoaz,aoast,       ...
                 xas,yas,Psi,igen,ipit,iyaw);

Res = Lpsi(dret,dret)\(Rpsi(dret) + Rgrav(dret));
Rval = sqrt((Res.')*Res);

dxpsi(dret) = Res;

%[dret Res.*Res]
%return
% 1:104    modal positions
% 105:208  modal velocities
% 209:250  collective aero
% 251:292  cos aero
% 293:334  sin aero
% 335:340  pitch
% 341:342  yaw
% 343:367  electrical
% 368:398  control

printf('   Initial residual %+5.3e\n',Rval);
fflush(stdout);

% ===================================================================
% Newton-Raphson method solution.
cnv = eps^0.3; % eps^0.6;
Ns = 100;
%bta = ones(Ns,1);
%bta = [0.01*ones(2,1);0.02*ones(2,1);0.05*ones(2,1);0.1*ones(2,1); ...
%       0.2*ones(2,1);0.5*ones(2,1);ones(Ns,1)];
%bta = [0.001*ones(2,1);0.002*ones(2,1);0.004*ones(2,1);0.008*ones(2,1); ...
%       0.016*ones(2,1);0.032*ones(2,1);0.064*ones(2,1);0.128*ones(2,1); ...
%       0.256*ones(2,1);0.5*ones(Ns,1)];
%bta = [0.01*ones(5,1);0.02*ones(5,1);0.04*ones(5,1);0.1*ones(5,1); ...
%       0.2*ones(5,1);0.5*ones(5,1);ones(Ns,1)];
bta = [0.01;0.02;0.05;0.1;0.2;0.5;ones(Ns,1)];
%bta = [0.01;0.02;0.05;0.1;0.2;0.5*ones(Ns,1)];
%bta = [0.0001;0.0002;0.0005;0.001;0.002;0.005;0.01;0.02;0.05;0.1;0.2;0.5;ones(Ns,1)];
%bta = 0.05*ones(Ns,1);
litmax = 10;
conv = 0;
iter = 0;
Niter = 0;
lam = 1;
while (((real(Rval) > cnv) && (iter < Ns)) || (iter == 0))

   iter = iter + 1;

   printf('   Iteration %5d, computing the tangent dynamics.\n',iter);
   fflush(stdout);

   % Compute the tangent function at the latest point.
   [Lpsi,Rpsi,ypsi,Apsi,Bpsi,Cpsi,Dpsi] =                   ...
            MBCCLT (1,xpsi,dxpsi,upsi,s,a,                  ...
                    epar,ppar,ypar,c,grav,P,shape0,mdamp0,  ...
                    Tas,Try,ch,Lel,foilwt,aoaz,aoast,       ...
                    xas,yas,Psi,igen,ipit,iyaw);
   dRdx = Lpsi(dret,dret)\Apsi(dret,dret);
%[slap,shp,ifrq] = eigVal (Lpsi(dret,dret)\Apsi(dret,dret));
%return
   printf('   Matrix condition %+6.3e\n',cond(dRdx));
   fflush(stdout);

   dxr = -dRdx\Res;
   lflg = 0;
   liter = 0;
   if (lam <= 0.5)
      lam = 2*lam;
   else
      lam = 1;
   end
   while ((lflg == 0) && (liter < litmax))
      liter = liter + 1;

      x1 = xpsi;
      x1(dret) = xpsi(dret) + bta(iter)*lam*dxr;

      % (Recall that dxpsi does not directly influence the nonlinear
      %  outputs, so I can use the existing value without updating.)
      [L1psi,R1psi,y1psi,jApsi,jBpsi,jCpsi,jDpsi] =            ...
               MBCCLT (0,x1,dxpsi,upsi,s,a,                    ...
                       epar,ppar,ypar,c,grav,P,shape0,mdamp0,  ...
                       Tas,Try,ch,Lel,foilwt,aoaz,aoast,       ...
                       xas,yas,Psi,igen,ipit,iyaw);
      Res1 = L1psi(dret,dret)\(R1psi(dret) + Rgrav(dret));
      R1 = sqrt((Res1.')*Res1);

      Rvec = Res1.*Res1;
      [Rmax,iRm] = max(sqrt(Rvec));

      ifndF = 8;
      iflp = imdofs(6) + 1;
      iedg = iflp + 1;

      xaero = Psi*TpsixB(2*Neta+[1:Nxa],:)*x1;
      printf('   %5d %5d  %+5.3e  %+5.3e\n',iter,liter,Rval,R1);
      printf('   Max residual %+5.3e, DOF %8d\n',Rmax,dret(iRm));
      printf('   W: %10.4f, b: %10.4f, Y: %10.4f\n', ...
             x1(2*Neta-4),x1(Neta-2),x1(Neta-5));
      printf('   qF %10.4f, qFd %10.4f, qf %10.4f, qfd %10.4f, qe %10.4f, qed %10.4f\n', ...
             x1(ifndF),x1(Neta+ifndF),x1(iflp),x1(Neta+iflp),x1(iedg),x1(Neta+iedg));
naoff = 10;
inda = 7*(naoff-1);
      printf('   ad %8.4f, a1 %8.4f, a2 %8.4f, Vihz %8.4f, Viht %8.4f, Viz %8.4f, Vit %8.4f\n', ...
             xaero(inda+1),xaero(inda+2),xaero(inda+3),xaero(inda+4), ...
             xaero(inda+5),xaero(inda+6),xaero(inda+7));
%      printf('   ad0 %10.4f, Viz0 %10.4f, adc %10.4f, Vizc %10.4f, ads %10.4f, Vizs %10.4f\n', ...
%             x1(2*Neta+7*naoff+1),x1(2*Neta+7*naoff+6), ...
%             x1(2*Neta+(Naero/3)+7*naoff+1),x1(2*Neta+(Naero/3)+7*naoff+6), ...
%             x1(2*Neta+2*(Naero/3)+7*naoff+1),x1(2*Neta+2*(Naero/3)+7*naoff+6));
      printf('   PehRSC %10.4f, Pem %10.4f, igq %10.4f, Vdc %+10.4f\n', ...
             x1(2*Neta+Naero+8+25+7),x1(2*Neta+Naero+8+25+15), ...
             x1(2*Neta+Naero+8+2),x1(2*Neta+Naero+8+3));
      fflush(stdout);

      if (real(R1) < Rval)
         % OK!  Prepare for the next iteration.
         lflg = 1;
         Res = Res1;
         Rval = R1;
         xpsi = x1;
         dxpsi(dret) = Res1;
      else
         % Backtrack.
         lam = 0.5*lam;
         if (liter == litmax)
            [iter R1 Rval]
            printf('Warning, proceeding without lambda convergence.\n');
            lflg = 1;
            Res = Res1;
            Rval = R1;
            xpsi = x1;
            dxpsi(dret) = Res1;
return
         end
      end

   end  % Newton inner.

   if (iter == Ns) && (real(Rval) > cnv)
      Rval
      printf('Warning, max iterations, proceeding without Rval convergence.\n');
return
   end

end

% Compute the final values.
dxpsi(dret) = Lpsi(dret,dret)\Rpsi(dret);
[Lpsi,Rpsi,ypsi,Apsi,Bpsi,Cpsi,Dpsi] =                   ...
         MBCCLT (1,xpsi,dxpsi,upsi,s,a,                  ...
                 epar,ppar,ypar,c,grav,P,shape0,mdamp0,  ...
                 Tas,Try,ch,Lel,foilwt,aoaz,aoast,       ...
                 xas,yas,Psi,igen,ipit,iyaw);

txt = outnm;

if (Vmag >= 10)
   Vstr = ['V' int2str(round(10*Vmag))];
else
   Vstr = ['V0' int2str(round(10*Vmag))];
end

eval(["save('-ascii','xpsi" txt Vstr ".txt','xpsi');"]);
eval(["save('-ascii','upsi" txt Vstr ".txt','upsi');"]);
eval(["save('-ascii','dxpsi" txt Vstr ".txt','dxpsi');"]);
eval(["save('-ascii','ypsi" txt Vstr ".txt','ypsi');"]);
eval(["save('-ascii','Rgrav" txt Vstr ".txt','Rgrav');"]);
eval(["save('-binary','Lpsi" txt Vstr ".bin','Lpsi');"]);
eval(["save('-binary','Apsi" txt Vstr ".bin','Apsi');"]);
eval(["save('-binary','Bpsi" txt Vstr ".bin','Bpsi');"]);
eval(["save('-binary','Cpsi" txt Vstr ".bin','Cpsi');"]);
eval(["save('-binary','Dpsi" txt Vstr ".bin','Dpsi');"]);
eval(["save('-binary','dret" txt Vstr ".bin','dret');"]);
eval(["save('-binary','shape" txt Vstr ".bin','shape0');"]);
eval(["save('-binary','mdamp" txt Vstr ".bin','mdamp0');"]);
eval(["save('-binary','bldof" txt Vstr ".bin','bldof');"]);

%geteigs;

