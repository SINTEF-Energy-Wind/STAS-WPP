% Run a Newton's method solution of the closed-loop turbine, and check the
% control tuning.  This attempts a single-level solution, that is, calling
% one function that returns dx/dt and A for the entire turbine.  It is also
% possible to pursue a hierarchical solution where each module is solved
% independently with inner Newton loops, and then the overall system with
% an outer Newton loop acting on the residuals in the interface variables.

clear;

printf('Reading inputs\n');
fflush(stdout);

nm = 'DTU10MW';

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

psiFlag = 1;   % = 1: MBC Vi formulation.
modFlag = 1;   % = 1: Apply body mode reduction to the structure.
shpFlag = 0;   % = 0: Use the input body mode basis functions for linearization.

% For linear analysis, the azimuth angle is arbitrary, as its
% influence will be eliminated by the MBC transform.
azi0 = 0*pi/180;
cp0 = cos(azi0);
sp0 = sin(azi0);

% ===================================================================
% Input parameters defining the load case.
a.dens = 1.225                     / ndens;
a.visc = 1.789e-5                  / nvisc;
Vmag   = 10                        / velocity;
yaw    = 0*pi/180;
alf    = 0.17;
grav   = [0;0;-9.807]              / (length/(time^2));
%bguess = max(1.5*(Vmag*velocity-10)*(pi/180),0);
bguess = 4*pi/180;

vs     = [33000;0]                 / voltage;
we     = 50*(2*pi)                 * time;
th_e   = 0;
Vhdc   = c.Vhdc;
Qh     = 0                         / power;

Area   = pi*(c.Ro^2);
%Pguess = minc(0.45*0.5*a.dens*Area*(Vmag^3),c.Pr);
%Pc     = c.Pr;
Pc     = 6.e6                       / power;
Pguess = Pc;
[Wguess,dW] = gains1 (Vmag,c.WVTab);

wg     = 0.5*c.np*Wguess;
ihg    = [1/current;-Pguess/(3000/voltage)];

Vinf   = [Vmag*cos(yaw);Vmag*sin(yaw);0];
t0     = 0;
betas  = [bguess;bguess;bguess];

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
% Fill out the global velocity vector based on the inputs.
Vg = zeros(3*Nae,1);
Vg(1:3:3*Nae-2) = Vinf(1);
Vg(2:3:3*Nae-1) = Vinf(2);
Vg(3:3:3*Nae)   = Vinf(3);

% Initial guess for induced velocity.
Vi0 = [-0.35*Vinf(1);0];
Viguess = zeros(2*Nae,1);
Viguess(1:2:2*Nae-1) = Vi0(1);
Viguess(2:2:2*Nae)   = Vi0(2);

% ===================================================================
% Initialize the structral calculation.
printf('Initializing structural states and mode shapes\n');
fflush(stdout);

[xs0,etas0,q0,dq0dt,d2q0dt2,P,shape0,shs,freq0,mdamp0,ret,slv, ...
 Ndj,Nnod,Neta] = structInit (s,yaw,azi0,betas,Wguess,modFlag);

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);
Nret = size(ret,1);
Nslv = size(slv,1);
Nylin = 4*Ndj + 9*Nnod + 137*Nae + 5;
Nu = Ndj + 3*Nae;

Ydof  = idofs(3);
Ddof  = idofs(4);
nodof = idofm(6) - 6;

iazi = Nylin - 4;
iW   = Nylin - 3;

% ===================================================================
% BEM setup, parameters that stay fixed throughout the calculations,
% and initial estimates of the aerodynamic states.
printf('Initializing aerodynamic states\n');
fflush(stdout);

[Tas,ch,Lel,foilwt,aoaz,aoast,xas,yas,iq] = BEMsetup (s,a);

Td_n = [cp0 -sp0 0;sp0 cp0 0;0 0 1];
Try = Tn_y;
Tyy0 = TFromTheta (q0(idofs(3)+[4:6]));
Ty0g = TFromTheta (P(idofs(3)+[4:6]));
Trg = Ty0g*Tyy0*Try;

[Tar,Trg,TB0g,TsB,TBB0,dTar,dTsB,dTBB0,wga] = ...
                      BEMprepTransforms (s,a,q0,dq0dt,P,Tas);

[zr,Area,Dp,rp,Lp,xeg,xhg,xyg] = ...
         BEMprepProjections (s,a,q0,dq0dt,P,Try,Trg);
Dps = zeros(Nae,1);
Dps(1:Neb)         = Dp(1);
Dps(Neb+[1:Neb])   = Dp(2);
Dps(2*Neb+[1:Neb]) = Dp(3);

xa0 = BEMinit (psiFlag,                          ...
               Viguess,Tar,Trg,Vg,wga,zr,ch,Lel, ...
               aoast,a.dens,Area,Dps,azi0,Wguess);

% Best estimate of the initial modal aero states.
bsh = bladeModeShape (s,ret,shape0);
Psi = aeroPsi (a,rp,bsh);
etaa0 = (Psi.')*xa0;
Naero = size(Psi,2);

% ===================================================================
% Actuators.
printf('Initializing actuator states\n');
fflush(stdout);

yp0 = zeros(9,1);
yp0([1 3 5]) = betas;
yp0([2 4 6]) = betas;
xp0 = zeros(6,1);
xp0([1 3 5]) = betas;
yy0 = zeros(3,1);
xy0 = zeros(2,1);

% ===================================================================
% Initialize the control states.
printf('Initializing control states\n');
fflush(stdout);

bavg = sum(betas)/3;

xc0 = [Wguess;Vinf(1);Wguess;bavg;bavg;Pguess;Pguess;Wguess;Wguess;bavg;bavg;bavg; ...
       0;0;Pguess;ihg(2);0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];


% ===================================================================
% Degrees-of-freedom to eliminate from the model.  In particular, 
% states associated with deactivated control functions must be
% eliminated.  Also the rotor azimuth should be eliminated.
%deldofs = [Neta-4, 2*Neta+Naero+8+25+[8 9 [12:24]]].';
%deldofs = Neta-4;

% ===================================================================
% Initial call to CLT function.  Get MBC dofs.
printf('Initial call to nonlinear turbine function\n');
fflush(stdout);

igen = [idofs(3)+[1:6] idofm(5)+[1:6] idofs(4)+[1:6] idofs(5)+[1:6]].';
ipit = [Ndof+[4:6]].';
iyaw = Ndof + 1;

u0 = [zeros(Ndj,1);Vg;we;th_e;vs;Pc;Qh];

% Initial states.
x0 = [etas0;etaa0;xp0;xy0;xe0;xc0];    % Guess

[Lcl,Rcl,yout,Acl,Bcl,Ccl,Dcl,blxdof,bludof,blydof] =                        ...
             buildClosedLoopTurbine (x0,u0,s,a,epar,ppar,ypar,c,             ...
                                     0,psiFlag,modFlag,shpFlag,              ...                                      
                                     grav,P,ret,slv,shape0,mdamp0,           ...
                                     Tas,Try,ch,Lel,foilwt,aoaz,aoast,       ...
                                     xas,yas,Psi,igen,ipit,iyaw);
Nx = size(Rcl,1);
Nu = size(Bcl,2);
Ny = size(yout,1);

% ===================================================================
% MBC transforms.
printf('Building MBC transforms.\n');
fflush(stdout);

[TpsixB,TBxpsi]   = MBC      (Nx,blxdof(:,1),blxdof(:,2),blxdof(:,3),azi0);
[TpsiuB,TBupsi]   = MBC      (Nu,bludof(:,1),bludof(:,2),bludof(:,3),azi0);
[TpsiyB,TBypsi]   = MBC      (Ny,blydof(:,1),blydof(:,2),blydof(:,3),azi0);
[dTpsixB,dTBxpsi] = derivMBC (Nx,blxdof(:,1),blxdof(:,2),blxdof(:,3),azi0);

%TpsixB = speye(Nx);
%TBxpsi = speye(Nx);
%TpsiuB = speye(Nu);
%TBupsi = speye(Nu);
%TpsiyB = speye(Ny);
%TBypsi = speye(Ny);
%dTpsixB = sparse(Nx,Nx);
%dTBxpsi = sparse(Nx,Nx);

trialCLT2;
