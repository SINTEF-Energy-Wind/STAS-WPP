clear;

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

% ====================================================================
% Trial calculation, decay test.  Calibrate full-model damping.  Get
% standstill frequencies (compare with linear model).
% 2 el: [8.17, 51.7] Hz.
runt = 0;
if (runt == 1)
   WW0 = 0;  % Overwrite, zero rotation.
   q0    = zeros(Ndof,1);
   dq0dt = zeros(Ndof,1);
   q0(i6c) = 0.02*(Pn(i6a).^2);
   x0 = [q0;dq0dt;azi0];
   ts = [0:0.00001:0.001].'; % [0:0.0005:1].';
   lsode_options('integration method','adams');
   [xs,ist,msg] = lsode (@(x,t)                                         ...
                         RotCantTimestep (x,t,mes,kes,PB,Pn,WW0,adamp), ...
                         x0,ts);
   Nt = size(ts,1);
   fid = fopen('out.txt','w');
   for it = 1:Nt
      fprintf(fid,'%+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e\n', ...
              ts(it),xs(it,1:3),xs(it,7:9),xs(it,Nx));
   end
   fclose('all');
end

% ====================================================================
% Task 1: Linear modal analysis, no rotation.  Compare with time-
% domain result.  2 el: [8.158, (32.631), 51.531] Hz.
run1 = 1;
if (run1 == 1)
   q0      = zeros(NdofB,1);
   dq0dt   = zeros(NdofB,1);
   d2q0dt2 = zeros(NdofB,1);
   P = [PB;Pn];
   [M,C,K] = RotCantLinFull (q0,dq0dt,d2q0dt2,P,mes,kes,adamp);
   A = [eye(Ndof) zeros(Ndof);zeros(Ndof) M] ...
     \ [zeros(Ndof) eye(Ndof);-K -C];
   [slap,shp,ifrq] = eigVal (A);
end

% ====================================================================
% Task 2: Nonlinear, full-model structural simulation in the time
% domain, with rotation, using direct equations.  Test nonlinear
% equations.  2 el, gam=2: [9.485,52.63] Hz.
tic
run2 = 1;
if (run2 == 1)
   q0    = zeros(Ndof,1);
   dq0dt = zeros(Ndof,1);  % Ramp up inside RotCantTimestep.
   x0 = [q0;dq0dt;azi0];
   ts = [0:0.001:0.1].'; 
   lsode_options('integration method','adams');
   [xs,ist,msg] = lsode (@(x,t)                                            ...
                         RotCantTimestep (x,t,mes,kes,PB,Pn,               ...
                                          feval('RCW',t,0.05,WW0),adamp,   ...
                                          feval('RCF',t,0,0.02,100,Ndof)), ...
                         x0,ts);
   Nt = size(ts,1);
   fid = fopen('outA.txt','w');
   for it = 1:Nt
      fprintf(fid,'%+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e\n', ...
              ts(it),xs(it,1:3),xs(it,7:9),xs(it,Nx));
   end
   fclose('all');
end
toc

% ====================================================================
% Task 3: Quasi-static solution of deformation under constant
% rotational speed, in the rotating frame.
run3 = 0;
if (run3 == 1)

   q0    = zeros(Ndof,1);
   dq0dt = zeros(Ndof,1);

   % Iterate to find the converged solution for the displacement.
   bta = 1;
   alphas = 1; % [0.5 1].';
   cnv = 1e-6;
   Na = size(alphas,1);
   Ns = 100;
   conv = zeros(Na,1);
   Niter = zeros(Na,1);
   xs = zeros(Nx,Na+1);
   xs(:,1) = [q0;dq0dt;azi0];
   for ia = 1:Na
      loadfrac = alphas(ia);
      WWf = WW0*loadfrac;
      iter = 0;
      Rval = 1;
      while ((Rval > cnv) && (iter <= Ns))
         
         iter = iter + 1;

         % Compute the RHS at the given x.  10.0: a late enough time
         % that the WW function in RotCantTimestep has ramped up to 
         % full speed.
         RHS = RotCantTimestep (xs(:,ia),10.0,mes,kes,PB,Pn,WWf,adamp);

         % Form the A matrix.
         P = [PB;Pn];
         qq = [qB;xs(1:Ndof,ia)];
         dqqdt = [dqBdt;xs(Ndof+[1:Ndof],ia)];
         d2qqdt2 = [zeros(6,1);RHS(Ndof+[1:Ndof])];
         [M,C,K] = RotCantLinFull (qq,dqqdt,d2qqdt2,P,mes,kes,adamp);

         A = [eye(Ndof) zeros(Ndof);zeros(Ndof) M] ...
           \ [zeros(Ndof) eye(Ndof);-K -C];

         % Perturbation to x.  (Fixed azi.)
         dx = -A\RHS(1:2*Ndof);

         % Update x.  (Fixed azi.)
         xs(:,ia) = xs(:,ia) + bta*[dx;0]; 

         % Residual.  (Not including azi.)
         Rval = max(abs(RHS(1:2*Ndof)));

      end

      if (Rval <= cnv)
         conv(ia) = 1;
      end
      Niter(ia) = iter;

      xs(:,ia+1) = xs(:,ia);

   end

end


% ====================================================================
% Task 4: Compare the linear and nonlinear solutions for a single
% element, given the initial stretch from rotation.  (Need run3 = 1.)
%
% Nonlinear equation of motion:
% M*d2q/dt2 = -(C+G)*dq/dt + H - K + Q*F
% Linear equation of motion:
% M*d2q/dt2 = -(dG/dqd qd0 + G0 - dH/dqd + dC/dqd qd0 + C0)*dq/dt
%     - (dM/dq qdd0 + dG/dq qd0 - dH/dq + dC/dq qd0 + dK/dq - dQ/dq F0)*q
%     + Q0 F
%     = - C' dq/dt - K' q + Q0 F 
run4 = 0;
if (run4 == 1)

   iel = 2;
   ic6 = 6*(iel-1);

   x0 = xs(:,Na);

   P = [PB;Pn];
   qq = [qB;x0(1:Ndof)];
   dqqdt = [dqBdt;x0(Ndof+[1:Ndof])];
   d2qqdt2 = [zeros(6,1);RHS(Ndof+[1:Ndof])];

   qn1      = qq(ic6+[1:6]);
   dqn1dt   = dqqdt(ic6+[1:6]);
   d2qn1dt2 = d2qqdt2(ic6+[1:6]);
   Pn1      = P(ic6+[1:6]);
   qn2      = qq(ic6+[7:12]);
   dqn2dt   = dqqdt(ic6+[7:12]);
   d2qn2dt2 = d2qqdt2(ic6+[7:12]);
   Pn2      = P(ic6+[7:12]);

   dqedt   = [dqBdt;dqn1dt;dqn2dt];
   d2qedt2 = [zeros(6,1);d2qn1dt2;d2qn2dt2];

   Qu1   = Qunod   (qn1,qB,Pn1,PB);
   Qu2   = Qunod   (qn2,qB,Pn2,PB);
   dQu1  = dQudq   (qn1,qB,Pn1,PB);
   dQu2  = dQudq   (qn2,qB,Pn2,PB);
   d2Qu1 = d2Qudq2 (qn1,qB,Pn1,PB);
   d2Qu2 = d2Qudq2 (qn2,qB,Pn2,PB);

   [xe,TsB] = elementCSFromNodes   (qn1,qn2,Pn1,Pn2);
   dTsB     = derivElementCS       (qn1,qn2,Pn1,Pn2,TsB);
   d2TsB    = secondDerivElementCS (qn1,qn2,Pn1,Pn2,TsB);

   mu   = getMu   (qn1,qn2,Pn1,Pn2,TsB);
   dmu  = dmudq   (qn1,qn2,Pn1,Pn2,TsB,dTsB);
   d2mu = d2mudq2 (qn1,qn2,Pn1,Pn2,TsB,dTsB,d2TsB);

   [mme,cce,kke] = buildElementMatrices (mes,kes,dqedt,d2qedt2,         ...
                                         Qu1,Qu2,dQu1,dQu2,d2Qu1,d2Qu2, ...
                                         mu,dmu,d2mu,TsB,dTsB,d2TsB);

   [cde,kde] = dampingCLin (adamp*kes,dqedt,dmu,d2mu);

   % Perturb the DOFs in turn.
   del = sqrt(eps);
   mmec = zeros(18,18);
   ccec = zeros(18,18);
   kkec = zeros(18,18);
   cdec = zeros(18,18);
   kdec = zeros(18,18);
   for idof = 1:6

      qBc = qB;
      qBc(idof) = qBc(idof) + i*del;
      dqedt   = [dqBdt;dqn1dt;dqn2dt];
      d2qedt2 = [zeros(6,1);d2qn1dt2;d2qn2dt2];
      Qu1c   = Qunod   (qn1,qBc,Pn1,PB);
      Qu2c   = Qunod   (qn2,qBc,Pn2,PB);
      dQu1c  = dQudq   (qn1,qBc,Pn1,PB);
      dQu2c  = dQudq   (qn2,qBc,Pn2,PB);
      [xe,TsB] = elementCSFromNodes   (qn1,qn2,Pn1,Pn2);
      dTsB     = derivElementCS       (qn1,qn2,Pn1,Pn2,TsB);
      mu   = getMu   (qn1,qn2,Pn1,Pn2,TsB);
      dmu  = dmudq   (qn1,qn2,Pn1,Pn2,TsB,dTsB);

      [mmc,ggc,hhc,kkc] = buildElementNL (mes,kes,dqedt,         ...
                                          Qu1c,Qu2c,dQu1c,dQu2c, ...
                                          mu,dmu,TsB,dTsB);

      cdc = dampingC (adamp*kes,dmu);
      kkec(:,idof) = (imag(mmc)/del)*d2qedt2 ...
                   + (imag(ggc)/del)*dqedt   ...
                   -  imag(hhc)/del          ...
                   +  imag(kkc)/del;

      qn1c = qn1;
      qn1c(idof) = qn1c(idof) + i*del;
      dqedt   = [dqBdt;dqn1dt;dqn2dt];
      d2qedt2 = [zeros(6,1);d2qn1dt2;d2qn2dt2];
      Qu1c   = Qunod   (qn1c,qB,Pn1,PB);
      Qu2    = Qunod   (qn2,qB,Pn2,PB);
      dQu1c  = dQudq   (qn1c,qB,Pn1,PB);
      dQu2   = dQudq   (qn2,qB,Pn2,PB);
      [xec,TsBc] = elementCSFromNodes   (qn1c,qn2,Pn1,Pn2);
      dTsBc     = derivElementCS       (qn1c,qn2,Pn1,Pn2,TsBc);
      muc   = getMu   (qn1c,qn2,Pn1,Pn2,TsBc);
      dmuc  = dmudq   (qn1c,qn2,Pn1,Pn2,TsBc,dTsBc);
      [mmc,ggc,hhc,kkc] = buildElementNL (mes,kes,dqedt,         ...
                                          Qu1c,Qu2,dQu1c,dQu2, ...
                                          muc,dmuc,TsBc,dTsBc);
      cdc = dampingC (adamp*kes,dmuc);
      kkec(:,idof+6) = (imag(mmc)/del)*d2qedt2 ...
                     + (imag(ggc)/del)*dqedt   ...
                     -  imag(hhc)/del          ...
                     +  imag(kkc)/del;

      qn2c = qn2;
      qn2c(idof) = qn2c(idof) + i*del;
      dqedt   = [dqBdt;dqn1dt;dqn2dt];
      d2qedt2 = [zeros(6,1);d2qn1dt2;d2qn2dt2];
      Qu1    = Qunod   (qn1,qB,Pn1,PB);
      Qu2c   = Qunod   (qn2c,qB,Pn2,PB);
      dQu1   = dQudq   (qn1,qB,Pn1,PB);
      dQu2c  = dQudq   (qn2c,qB,Pn2,PB);
      [xec,TsBc] = elementCSFromNodes   (qn1,qn2c,Pn1,Pn2);
      dTsBc      = derivElementCS       (qn1,qn2c,Pn1,Pn2,TsBc);
      muc   = getMu   (qn1,qn2c,Pn1,Pn2,TsBc);
      dmuc  = dmudq   (qn1,qn2c,Pn1,Pn2,TsBc,dTsBc);
      [mmc,ggc,hhc,kkc] = buildElementNL (mes,kes,dqedt,       ...
                                          Qu1,Qu2c,dQu1,dQu2c, ...
                                          muc,dmuc,TsBc,dTsBc);
      cdc = dampingC (adamp*kes,dmuc);
      kkec(:,idof+12) = (imag(mmc)/del)*d2qedt2 ...
                      + (imag(ggc)/del)*dqedt   ...
                      -  imag(hhc)/del          ...
                      +  imag(kkc)/del;

      dqBdtc = dqBdt;
      dqBdtc(idof) = dqBdtc(idof) + i*del;
      dqedtc   = [dqBdtc;dqn1dt;dqn2dt];
      d2qedt2 = [zeros(6,1);d2qn1dt2;d2qn2dt2];
      Qu1   = Qunod   (qn1,qB,Pn1,PB);
      Qu2   = Qunod   (qn2,qB,Pn2,PB);
      dQu1  = dQudq   (qn1,qB,Pn1,PB);
      dQu2  = dQudq   (qn2,qB,Pn2,PB);
      [xe,TsB] = elementCSFromNodes   (qn1,qn2,Pn1,Pn2);
      dTsB     = derivElementCS       (qn1,qn2,Pn1,Pn2,TsB);
      mu   = getMu   (qn1,qn2,Pn1,Pn2,TsB);
      dmu  = dmudq   (qn1,qn2,Pn1,Pn2,TsB,dTsB);

      [mmc,ggc,hhc,kkc] = buildElementNL (mes,kes,dqedtc,    ...
                                          Qu1,Qu2,dQu1,dQu2, ...
                                          mu,dmu,TsB,dTsB);

      cdc = dampingC (adamp*kes,dmu);
      ccec(:,idof) =  real(ggc(:,idof))            ...
                   + (imag(ggc)/del)*real(dqedtc)  ...
                   -  imag(hhc)/del;

      dqn1dtc = dqn1dt;
      dqn1dtc(idof) = dqn1dtc(idof) + i*del;
      dqedtc   = [dqBdt;dqn1dtc;dqn2dt];
      d2qedt2 = [zeros(6,1);d2qn1dt2;d2qn2dt2];
      Qu1    = Qunod   (qn1,qB,Pn1,PB);
      Qu2    = Qunod   (qn2,qB,Pn2,PB);
      dQu1   = dQudq   (qn1,qB,Pn1,PB);
      dQu2   = dQudq   (qn2,qB,Pn2,PB);
      [xe,TsB] = elementCSFromNodes   (qn1,qn2,Pn1,Pn2);
      dTsB     = derivElementCS       (qn1,qn2,Pn1,Pn2,TsB);
      mu   = getMu   (qn1,qn2,Pn1,Pn2,TsB);
      dmu  = dmudq   (qn1,qn2,Pn1,Pn2,TsB,dTsB);
      [mmc,ggc,hhc,kkc] = buildElementNL (mes,kes,dqedtc,    ...
                                          Qu1,Qu2,dQu1,dQu2, ...
                                          mu,dmu,TsB,dTsB);
      cdc = dampingC (adamp*kes,dmu);
      ccec(:,idof+6) =  real(ggc(:,idof+6))        ...
                   + (imag(ggc)/del)*real(dqedtc)  ...
                   -  imag(hhc)/del;

      dqn2dtc = dqn2dt;
      dqn2dtc(idof) = dqn2dtc(idof) + i*del;
      dqedtc   = [dqBdt;dqn1dt;dqn2dtc];
      d2qedt2 = [zeros(6,1);d2qn1dt2;d2qn2dt2];
      Qu1    = Qunod   (qn1,qB,Pn1,PB);
      Qu2    = Qunod   (qn2,qB,Pn2,PB);
      dQu1   = dQudq   (qn1,qB,Pn1,PB);
      dQu2   = dQudq   (qn2,qB,Pn2,PB);
      [xe,TsB] = elementCSFromNodes   (qn1,qn2,Pn1,Pn2);
      dTsB     = derivElementCS       (qn1,qn2,Pn1,Pn2,TsB);
      mu   = getMu   (qn1,qn2,Pn1,Pn2,TsB);
      dmu  = dmudq   (qn1,qn2,Pn1,Pn2,TsB,dTsB);
      [mmc,ggc,hhc,kkc] = buildElementNL (mes,kes,dqedtc,    ...
                                          Qu1,Qu2,dQu1,dQu2, ...
                                          mu,dmu,TsB,dTsB);
      cdc = dampingC (adamp*kes,dmu);
      ccec(:,idof+12) =  real(ggc(:,idof+12))      ...
                   + (imag(ggc)/del)*real(dqedtc)  ...
                   -  imag(hhc)/del;

   end

%full(kke-kkec)
full(cce-ccec)

end

% ====================================================================
% Task 5: Linear modal analysis, rotating.  Initial condition from
% quasi-static solution.  Test linear equations.  (Need run3 = 1.)
% 2 el, gam=2: [9.585,51.87] Hz.
run5 = 0;
if (run5 == 1)

   x0 = xs(:,Na);

   P = [PB;Pn];
   qq = [qB;x0(1:Ndof)];
   dqqdt = [dqBdt;x0(Ndof+[1:Ndof])];
   d2qqdt2 = [zeros(6,1);RHS(Ndof+[1:Ndof])];
   [M,C,K] = RotCantLinFull (qq,dqqdt,d2qqdt2,P,mes,kes,adamp);
   A = [eye(Ndof) zeros(Ndof);zeros(Ndof) M] ...
     \ [zeros(Ndof) eye(Ndof);-K -C];
   [slap,shp,ifrq] = eigVal (A);

   [[1:size(slap,1)].' imag(slap)*T]

end


 





