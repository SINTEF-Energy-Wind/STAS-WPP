% ===================================================================
% Initialize the MBC variables, and get the N-R residual.
printf('Newton-Raphson solution.\n');
printf('   Computing the initial residual.\n');
fflush(stdout);

%xpsi = TBxpsi*x0;                  % Based on initial guess
lcnm = 'xpsi_V100_P060_reduced';             % ... or read from file.
eval(['load "' lcnm '.txt";']);
eval(['xpsi = ' lcnm ';']);
%load 'xpsi.txt';

%xpsi0 = xpsi;
%xpsi = zeros(Nx,1);
%xpsi(1:2*Neta+Naero+8+25+11) = xpsi0(1:2*Neta+Naero+8+25+11);
%xpsi(2*Neta+Naero+8+25+12) = xpsi0(2*Neta+Naero+8+25+4);
%xpsi(2*Neta+Naero+8+25+[13:31]) = xpsi0(2*Neta+Naero+8+25+[12:30]);

%xpsi0 = xpsi;
%xpsi = zeros(Nx,1);
%xpsi([1:Neta]) = xpsi0([[1:6] 7 8 17 18 27 28 29 30 31 32 47 48 63 64 79:84]);
%xpsi(Neta+[1:Neta]) = xpsi0(84+[[1:6] 7 8 17 18 27 28 29 30 31 32 47 48 63 64 79:84]);
%xpsi(2*Neta+[1:7]) = xpsi0(168+7*4+[1:7]);
%xpsi(2*Neta+7+[1:7]) = xpsi0(168+42+7*4+[1:7]);
%xpsi(2*Neta+14+[1:7]) = xpsi0(168+84+7*4+[1:7]);
%xpsi(2*Neta+21+[1:8+25+31]) = xpsi0(168+126+[1:8+25+31]);

upsi = TBupsi*u0;

%xc = xpsi(2*Neta+Naero+8+25+[1:30]);
%xpsi(2*Neta+Naero+8+25+[1:30]) = ...
%       [xc(1);xc(2);xc(3);0.05;0.05;Pguess;Pguess;xc(8);xc(9);0.05;0.05; ...
%        xc(12);xc(13);Pguess;ihg(2);xc(16:30)];


%[jnkvec,dret,dslv] = partitionMatrix ([1:size(xpsi,1)].',deldofs,[]);

% Locked yaw, azimuth.  Only RSC, Vobs, gen P control.
%dret = [[1:Neta-6] Neta-[3 2 1 0] Neta+[1:Neta-6] 2*Neta-[4 3 2 1 0] ...
%        2*Neta+[1:Naero] 2*Neta+Naero+[1:8] 2*Neta+Naero+8+[1:25] ...
%        2*Neta+Naero+8+25+[1 2 3 4 5 6 7 8 9 10 11 12 15 16 17 18]].';

% Locked rotor speed.
%dret = [[1:Neta-6] Neta-[3 2 1 0] Neta+[1:Neta-6] 2*Neta-[3 2 1 0] ...
%        2*Neta+[1:Naero] 2*Neta+Naero+[1:8] 2*Neta+Naero+8+[1:25] ...
%        2*Neta+Naero+8+25+[1 2 3 4 5 6 7 8 9 10 11 12 15 16 17 18]].';

% Rigid structure.
%dret = [Neta-[3 2 1 0] 2*Neta-[4 3 2 1 0] ...
%        2*Neta+[1:Naero] 2*Neta+Naero+[1:8] 2*Neta+Naero+8+[1:25] ...
%        2*Neta+Naero+8+25+[1 2 3 4 5 6 7 8 9 10 11 12 15 16 17 18]].';

% Collective only.
dret = [[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16] Neta-[3 2 1 0] ...
        Neta+[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16] 2*Neta-[4 3 2 1 0] ...
        2*Neta+[1:Naero/3] 2*Neta+Naero+[1:2] 2*Neta+Naero+8+[1:25] ...
        2*Neta+Naero+8+25+[1 2 3 4 5 6 7 8 9 10 11 12 15 16 17 18]].';

c = STASControl_DTU10MW ();

[Lpsi,Rpsi,ypsi,Apsi,Bpsi,Cpsi,Dpsi] =                ...
         MBCCLT (xpsi,upsi,Neta,                      ...
                 TpsixB,dTpsixB,TBxpsi,TpsiuB,TBypsi, ...
                 s,a,epar,ppar,ypar,c,                ...
                 0,psiFlag,modFlag,shpFlag,           ...                                      
                 grav,P,ret,slv,shape0,mdamp0,        ...
                 Tas,Try,ch,Lel,foilwt,aoaz,aoast,    ...
                 xas,yas,Psi,igen,ipit,iyaw);
Res = Lpsi(dret,dret)\Rpsi(dret);
Rval = (Res.')*Res;

%[dret Res.*Res]
%return
% 1:84     modal positions
% 85:168   modal velocities
% 169:210  collective aero
% 211:252  cos aero
% 253:294  sin aero
% 295:300  pitch
% 301:302  yaw
% 303:327  electrical
% 328:358  control

printf('   Initial residual %+5.3e\n',Rval);
fflush(stdout);

% ===================================================================
% Newton-Raphson method solution.
cnv = eps^0.2; % eps^0.6;
Ns = 100;
%bta = [0.01*ones(2,1);0.02*ones(2,1);0.04*ones(4,1);0.1*ones(5,1); ...
%       0.2;0.4;0.5*ones(4,1);ones(Ns,1)];
%bta = [0.001*ones(2,1);0.002*ones(2,1);0.004*ones(2,1);0.008*ones(2,1); ...
%       0.016*ones(2,1);0.032*ones(2,1);0.064*ones(2,1);0.128*ones(2,1); ...
%       0.256*ones(2,1);0.5*ones(Ns,1)];
%bta = [0.01*ones(5,1);0.02*ones(5,1);0.04*ones(5,1);0.1*ones(5,1); ...
%       0.2*ones(5,1);0.5*ones(10,1)];
%bta = [0.01;0.02;0.05;0.1;0.2;0.5;ones(Ns,1)];
bta = ones(Ns,1);
litmax = 20;
conv = 0;
iter = 0;
Niter = 0;
lam = 1;
while (((real(Rval) > cnv) && (iter < Ns)) || (iter == 0))
   iter = iter + 1;

   printf('   Iteration %5d, computing the tangent dynamics.\n',iter);
   fflush(stdout);

   % Compute the tangent function at the latest point.
   [Lpsi,Rpsi,ypsi,Apsi,Bpsi,Cpsi,Dpsi] =                  ...
            MBCCLT (xpsi,upsi,Neta,                      ...
                    TpsixB,dTpsixB,TBxpsi,TpsiuB,TBypsi, ...
                    s,a,epar,ppar,ypar,c,                ...
                    1,psiFlag,modFlag,shpFlag,           ...                                      
                    grav,P,ret,slv,shape0,mdamp0,        ...
                    Tas,Try,ch,Lel,foilwt,aoaz,aoast,    ...
                    xas,yas,Psi,igen,ipit,iyaw);
   dRdx = Lpsi(dret,dret)\Apsi(dret,dret);

   printf('   Matrix condition %+6.3e\n',cond(dRdx));
   fflush(stdout);

%[uu,ss,vv] = svd(dRdx);
%N = size(dRdx,1);
%[ss(N,N) uu(dret==332,N) vv(dret==332,N)]

%[slap,shp,ifrq] = eigVal (Apsi(dret,dret));

%return

   dxr = -dRdx\Res;
   lflg = 0;
   liter = 0;
   while ((lflg == 0) && (liter < litmax))
      liter = liter + 1;

      x1 = xpsi;
      x1(dret) = xpsi(dret) + bta(iter)*lam*dxr;

      [L1psi,R1psi,y1psi,jnkApsi,jnkBpsi,jnkCpsi,jnkDpsi] = ...
               MBCCLT (x1,upsi,Neta,                        ...
                       TpsixB,dTpsixB,TBxpsi,TpsiuB,TBypsi, ...
                       s,a,epar,ppar,ypar,c,                ...
                       0,psiFlag,modFlag,shpFlag,           ...                                      
                       grav,P,ret,slv,shape0,mdamp0,        ...
                       Tas,Try,ch,Lel,foilwt,aoaz,aoast,    ...
                       xas,yas,Psi,igen,ipit,iyaw);
      Res1 = L1psi(dret,dret)\R1psi(dret);
      R1 = (Res1.')*Res1;

      Rvec = Res1.*Res1;
      [Rmax,iRm] = max(Rvec);

      printf('   %5d %5d  %+5.3e  %+5.3e\n',iter,liter,Rval,R1);
      printf('   Max residual %+5.3e, DOF %8d\n',Rmax,dret(iRm));
      printf('   W: %10.4f, b: %10.4f\n', ...
             x1(2*Neta-4),x1(Neta-2));
%      printf('   qF %10.4f, qFd %10.4f, qf %10.4f, qfd %10.4f, qe %10.4f, qed %10.4f\n', ...
%             x1(8),x1(Neta+8),x1(31),x1(Neta+31),x1(32),x1(Neta+32));
      printf('   qF %10.4f, qFd %10.4f, qf %10.4f, qfd %10.4f, qe %10.4f, qed %10.4f\n', ...
             x1(8),x1(Neta+8),x1(15),x1(Neta+15),x1(16),x1(Neta+16));
naoff = 0;
      printf('   ad0 %10.4f, Viz0 %10.4f, adc %10.4f, Vizc %10.4f, ads %10.4f, Vizs %10.4f\n', ...
             x1(2*Neta+7*naoff+1),x1(2*Neta+7*naoff+6), ...
             x1(2*Neta+(Naero/3)+7*naoff+1),x1(2*Neta+(Naero/3)+7*naoff+6), ...
             x1(2*Neta+2*(Naero/3)+7*naoff+1),x1(2*Neta+2*(Naero/3)+7*naoff+6));
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

save('-ascii','xpsi.txt','xpsi');

%geteigs;
