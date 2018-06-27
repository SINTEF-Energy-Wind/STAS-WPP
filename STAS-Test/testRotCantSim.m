% Simulate a rotating cantilever with an implicit integration method.

clear;

tic

% Refer to Yoo and Shin (1998) for a definition of input parameters.
% T   = sqrt(rho*A*L^4 / EI)
% tau = t/T
% xi  = x/L
% th  = q/L   q: modal amplitude.
% del = r/L   r: root radius.
% gam = W*T
% wn  = w*T
% alf = sqrt(A*L^2 / Izz)

% Beam properties and operating condition.
del = 0; % 5;
gam = 0; % 10;
alf = 70;
Nel = 2;

Len = 1;
wid = sqrt(12)*Len/alf;
thk = wid*1.000001;  % For simplicity, but not precisely symmetric.
rho = 7850;
EE = 2e11;
GG = EE/2.6;
adamp = 0.00001; % 0.000005; %

Area = wid*thk;
IIyy = (wid*(thk^3))/12;
IIzz = (thk*(wid^3))/12;
JJ = ((thk^3)*wid/3)*(1 - (192*thk/((pi^5)*wid)) ...
   * sum((1./((2*[1:8] - 1).^5)).*tanh(pi*(2*[1:8] - 1)*wid/(2*thk))));
EIyy = EE*IIyy;
EIzz = EE*IIzz;
GJ = GG*JJ;
EA = EE*Area;
r = del*Len;
T = sqrt(rho*Area*(Len^4)/EIyy);
rhoA = rho*Area;
rhoJ = rho*JJ;

Lel = Len/Nel;

WW0 = gam/T;

Ndof  = 6*Nel;
NdofB = 6*Nel + 6;
Nx = 2*Ndof + 1;  % q, dq/dt, azi.

% Define the basic element mass and stiffness matrices.  This is the
% same for each element.
EEs = diag([EA;0;0;GJ;EIyy;EIzz]);
rhos = diag([rhoA;rhoA;rhoA;rhoJ;0;0]);

kes = buildkes (EEs,Lel);
mes = buildmes (rhos,Lel);

% Some convenient indices.
i6a = [1:6:6*Nel-5].';
i6b = [2:6:6*Nel-4].';
i6c = [3:6:6*Nel-3].';
i6d = [4:6:6*Nel-2].';
i6e = [5:6:6*Nel-1].';
i6f = [6:6:6*Nel].';

% Undeformed, unrotated position.
azi0 = 0;
ca0 = cos(azi0);
sa0 = sin(azi0);

TB0g = [ca0 -sa0 0; sa0 ca0 0; 0 0 1];
PB = zeros(6,1);
PB(4:6) = thetaFromT (TB0g);
PB(1:3) = TB0g*[r;0;0];

Pn = zeros(6*Nel,1);
Pn(i6a) = ([1:Nel].')*Lel;

% Get the steady-state velocity of the reference node generalized
% coordinates.
TBg = [ca0 -sa0 0; sa0 ca0 0; 0 0 1];

qB = zeros(6,1);
qB(4:6) = thetaFromT (TBg);
qB(1:3) = TBg*[r;0;0];

dTBg = [-sa0 -ca0 0; ca0 -sa0 0; 0 0 0];
Gz = TBg*(dTBg.');
Gx = zeros(3,3);
Gy = zeros(3,3);

[jnk,dT] = dTdth (qB(4:6));
Fx = dT(:,1:3)*(TBg.');
Fy = dT(:,4:6)*(TBg.');
Fz = dT(:,7:9)*(TBg.');

GG = [-Gx(2,3) -Gy(2,3) -Gz(2,3); ...
       Gx(1,3)  Gy(1,3)  Gz(1,3); ...
      -Gx(1,2) -Gy(1,2) -Gz(1,2)];
FF = [-Fx(2,3) -Fy(2,3) -Fz(2,3); ...
       Fx(1,3)  Fy(1,3)  Fz(1,3); ...
      -Fx(1,2) -Fy(1,2) -Fz(1,2)];

dqBdt = zeros(6,1);
dqBdt(4:6) = -FF\GG*[0;0;WW0];
dqBdt(1:3) = dTBg*[r;0;0]*WW0;

P = [PB;Pn];
q0    = zeros(Ndof,1);
dq0dt = zeros(Ndof,1);  % Ramp up inside RotCantTimestep.
x0 = [q0;dq0dt;azi0];
dt = 0.002;
ts = [0:dt:0.5].'; 
Nt = size(ts,1);

% Get an initial function value.
xn = x0;
xnp = xn;

xns = zeros(Nx,Nt);
Fs  = zeros(Ndof,Nt);

cnv   = eps^(0.6);
Ns    = 10;
bta   = [1.0*ones(1,100)].';
for it = 1:Nt

'================================='
ts(it)
   t = ts(it);

   % ===================================================================
   % Define functions of time here...
   t0 = 0;
   t1 = 0.1;
   Fmax = 100;
   F = RCF (t,t0,t1,Fmax,Ndof);
   Fs(:,it) = F;
   WW = RCW (t,0.05,WW0);  % Linear ramp up, > period of first mode.   
   % ===================================================================

   % Update the residual with the best guess, which is the function
   % evaluated at xn.  In the trapezoid case, this works out to
   % dxndt.
   dxndt = RotCantTimestep (xn,ts(1),mes,kes,PB,Pn,WW,adamp,F);
   dxnpdt = dxndt;
   Res = -dxndt;
   Rval = (Res.')*Res;
Rval
   conv  = 0;
   iter  = 0;
   Niter = 0;
   while ((Rval > cnv) && (iter < Ns))
      iter = iter + 1;

iter
Rval

      % Compute the tangent dynamics at the latest point.
      qq = [qB;xnp(1:Ndof)];
      dqqdt = [dqBdt;xnp(Ndof+[1:Ndof])];
      d2qqdt2 = [zeros(6,1);dxnpdt(Ndof+[1:Ndof])];
      [M,C,K] = RotCantLinFull (qq,dqqdt,d2qqdt2,P,mes,kes,adamp);
      A = [eye(Ndof) zeros(Ndof) zeros(Ndof,1); ...
           zeros(Ndof) M zeros(Ndof,1); ...
           zeros(1,2*Ndof) 1] ...
        \ [zeros(Ndof) eye(Ndof) zeros(Ndof,1); ...
           -K -C zeros(Ndof,1); ...
           zeros(1,2*Ndof+1)];
      dRdxn = speye(Nx)/dt - 0.5*A;

%full(dRdxn)
%Res

      % ----------------------------------------------------------------
      % Apply Newton's method on the MBC-transformed equations.
      lam = 1;
      dx = -dRdxn\Res;
      lflg = 0;
      litmax = 10;
      liter = 0;
      while ((lflg == 0) && (liter < litmax))
         liter = liter + 1;
'-----'
liter
         xnl = xnp + bta(iter)*lam*dx;
         dxnldt = RotCantTimestep (xnl,t,mes,kes,PB,Pn,WW,adamp,F);
         Res = (xnl - xn)/dt - 0.5*(dxnldt + dxndt);
         R1 = (Res.')*Res;
%[dx xn xnp xnl]
%[dxndt dxnpdt dxnldt Res]
[R1 Rval]
         if (R1 < Rval)
            % OK!  Prepare for the next iteration.
            lflg = 1;
            Rval = R1;
            xnp = xnl;
            dxnpdt = dxnldt;
         else
            % Backtrack.
            lam = 0.5*lam;
            if (liter == litmax)
               iter
               'Warning, proceeding without lambda convergence.'
               lflg = 1;
               xnp = xnl;
               dxnpdt = dxnldt;
            end
         end

      end  % Newton inner.

   end  % Newton outer.

   % Increment the n index in preparation for the next timestep.
   xn = xnp;
   xns(:,it) = xn;

%fidt = fopen ('stat.txt','a');
%fprintf(fidt,'%+5.6e %4d %+5.6e %+5.6e %+5.6e\n', ...
%        t,iter,Rval,F(Ndof-6+2),xn(Ndof-6+2));
%fclose('all');

end  % Timestepping.

fid = fopen ('outN.txt','w');
for it = 1:Nt
   fprintf(fid,'%+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e\n', ...
           ts(it),Fs(Ndof-6+2,it),xns(1:3,it),xns(7:9,it),xns(Nx,it));
end
fclose('all');

toc


