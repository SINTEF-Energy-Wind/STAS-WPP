% 05.08.2017  Complex step dS verified against finite difference.

clear;

[Nef,Net,Nev,Ner,Nen,Ned,Neh,Neb,icp,       ...
 Nmud,Nwater,Ntow,Nnac,Ndrv,Nbld,           ...
 Lf,Lt,Lv,Lr,Ln,Ld,Lh,Lb,                   ...
 Lel,xis,rhos,EEs,                          ...
 EE,GG,density,viscosity,                   ...
 wnod,rhow,CmA,Cdf0,Cdt0,                   ...
 Diaf,Diat,Diav,Diar,Dian,Diad,             ...
 xia,chord,xpc,foil,aoaz,z0,delta,phi,      ...
 kxg,kyg,kzg,kthzg,cxg,cyg,czg,cthzg,       ...
 Wsched,betaSched,Prated,Trated] = STASTurbine ();

load('-ascii','LTMnorms.txt');
length = LTMnorms(1);
time   = LTMnorms(2);
mass   = LTMnorms(3);
velocity = length/time;

epsilon = sqrt(eps);

flagg = 2;  % 0 for no gradient, 
            % 1 for gradient wrt R, 
            % 2 for gradient wrt Omega.

NN = 3*Neb;
Nel = Nef + Net + Nev + Ner + Nen + Ned + Neh + 3 + 3*Neb;

[nf,df,LCdef] = STASSetup ();
Nlc = size(LCdef,1);

towerFlag = 1;  % = 1: Include nodes on the tower in the turbulence
                %      cross-spectra.

psi0 = 0;
chi0 = 0;

Wr = Wsched(size(Wsched,1));

% Dummy call to get opFlag.
cpar = STASControl((10/velocity),Wsched,Wr,0,Prated);
opFlag = cpar(1);

Vinfs = LCdef(:,1);
Vmkss = Vinfs*velocity;
Pcs   = LCdef(:,2);
Vdc0s = LCdef(:,3);
vpd0s = LCdef(:,4);
we0s  = LCdef(:,5);
Is    = 0.14*(0.75*Vmkss + 5.6)./Vmkss;  % IEC 61400-1 NTM.
                                         % Class B, Iref = 0.14.
Lus   = (180/length)*ones(Nlc,1); % Based on surface roughness length of
                                  % about 0.01 m.  0.0001: 235, 0.001: 200,
                                  % 0.01: 175, 0.1: 150.  Intended to
                                  % represent a turbine deep in a plant.

z0 = 0.01/length;        % Surface roughness length for wind shear.

Ns = 3; % 8;                  % For "normal" wind shear, don't need more than 8 terms, 
                         % could be even fewer.  For wake-shadow, could use more.
Nt = 3; % 10;                 % Tower shadow.  Want to preserve maybe 30 terms, but
                         % 3*Nt*P should not exceed 2*nf*df.


Lbel = Lel(Nel-3*Neb+[1:Neb]');

r = zeros(Neb,1);
r(1) = Lh + 0.5*Lbel(1);
for iel = 2:Neb
   r(iel) = r(iel-1) + 0.5*(Lbel(iel) + Lbel(iel-1));
end
r3 = [r;r;r];

Ro = r(Neb) + 0.5*Lbel(Neb);

if (flagg == 1)
   Ro = Ro + i*epsilon;
   r = (r/real(Ro))*Ro;
   r3 = [r;r;r];
   Lh = (Lh/real(Ro))*Ro;
   Lb = Ro - Lh;
   ch = 'R';
end

Dia = 2*Ro;



if (towerFlag == 1)

   [Ogy_g,Oyd_y,Odp_d,Ogs_g,Oys_y,Ods_d,Ops_p] =     ...
      nodalOrigins (Nef,Net,Nev,Ner,Nen,Ned,Neh,Neb, ...
                    delta,Lf,Lt,Lv,Lr,Ln,Ld,Lh,Lb,Lel);

   [Tn_y,Th_d,Tb_h] = buildBasicTransforms (delta,phi);

   zet_g = zeros(Net,1);
   zn1 = 0;
   zn2 = Lel(Nef+1);
   for iel = 1:Net
      zet_g(iel) = 0.5*(zn1 + zn2);
      zn1 = zn1 + Lel(Nef+iel);
      zn2 = zn2 + Lel(Nef+iel+1);
   end

end

for ilc = 1:Nlc

   Vinf = Vinfs(ilc);
   Pc   = Pcs(ilc);
   Vdc0 = Vdc0s(ilc);
   vpd0 = vpd0s(ilc);
   we0  = we0s(ilc);
   I    = Is(ilc);
   Lu   = Lus(ilc);

   Vmks = Vinf*velocity;

   % Find the operating rotor speed Wop.
   Vlo = floor(Vmks);
   Vhi = ceil(Vmks);
   if (Vlo == Vhi)
      ilo = Vlo - 3;
      Wnom = Wsched(ilo);
      bnom = betaSched(ilo);
   else
      f = Vmks - Vlo;
      ilo = Vlo - 3;
      ihi = Vhi - 3;
      Wnom = (1-f)*Wsched(ilo) + f*Wsched(ihi);
      bnom = (1-f)*betaSched(ilo) + f*betaSched(ihi);
   end

   [Pnom,Pop,Wop,bop,Fapop,Virop] =                ...
             steadyOp (opFlag,Pc,Prated,Trated,    ...
                       Neb,Vinf,Wnom,bnom,Wr,      ...
                       psi0,Dia,density,viscosity, ...
                       Vdc0,vpd0,we0,              ...
                       xia,r,Lbel,chord,xpc,foil);

   ify = [2:6:6*NN-4]';
   ifz = [3:6:6*NN-3]';
   Fz = -Fapop(ify)*sin(bop(1)) + Fapop(ifz)*cos(bop(1));
   Ft = Fapop(ify)*cos(bop(1)) + Fapop(ifz)*sin(bop(1));

   Thrust = sum(Fz);
   Torque = sum(r3.*Ft);
   Pmech = Wop*Torque;
   beta0 = [bop;bop;bop];

real([Vinf Pnom Pop Pmech Torque Thrust Wop bop])

   if (flagg == 2)
      % Directly scaling W does not work, because this doesn't
      % correctly represent the "rolling" of the spectral peaks
      % along the frequency axis, as W changes.  Rather, scale
      % the time variable such that it follows W.  This implies
      % that W stays the same and the input Vinf gets scaled.
      % But only the convection velocity.  The amplitude of the
      % turbulent eddies should stay the same, so increase I
      % accordingly.
%      Wop = Wop + i*epsilon;
%      ch = 'W';
      Vinf = Vinf + i*epsilon;
      I = I*real(Vinf)/Vinf;
      ch = 'V';
   end

   [SpsiR,SpsiI] = mbcSij (Neb,nf,df,Vinf,I,Lu,Wop,r3,psi0);

   Vzrp = Vinf*ones(Neb,1);

   [SperR,SperI] =                                                   ...
          periodicWind (Nef,Net,Nev,Ner,Nen,Ned,Neh,Neb,Nmud,Nwater, ...
                        Lf,Lt,Lv,Lr,Ln,Ld,Lh,Lb,                     ...
                        Ns,Nt,Lel,Diaf,Diat,delta,phi,df,            ...
                        Wop,z0,Vinf,Vzrp);

   % The periodic wind terms contribute to 3nP multiples of the
   % rotational frequency.  It's a matter of finding the right
   % rows and columns in the Spsi matrix.
   frot = 3*(Wop/(2*pi))*[[-size(SperR,1):-1] [1:size(SperR,1)]].';
   f = df*[-2*nf+1:2*nf].';
   ir = interp1(f,[[2*nf+1:4*nf] [1:2*nf]].',real(frot),'nearest');
   ic = [1:4:4*((3*Neb)^2) - 3];

   % Eliminate phase information in the periodic signal.  In general, we
   % know its strength, but its phase would depend on the starting 
   % azimuth, which is not permissible for stochastic analysis.  Rather,
   % we preserve the strength of the signal and neglect the relative phase.
   SpsiR(ir,ic) = SpsiR(ir,ic) ...
                + sqrt(([SperR(size(SperR,1):-1:1,:);SperR]).^2 ...
                +      ([SperI(size(SperI,1):-1:1,:);SperI]).^2);

   clear SperR SperI;

   % Now compute and save the nominal spectra.
   Spsi    = zeros(size(SpsiR),'single');
   Spsi(:) = real(SpsiR);
   Spsi    = Spsi + i*real(SpsiI);

   % Memory management: we no longer need the real parts of SpsiR and
   % SpsiI.  Store the imaginary parts in SpsiI.
   SpsiI = imag(SpsiR) + i*imag(SpsiI);
   clear SpsiR;

   % Float halves the file size, and a comparison indicates no meaningful
   % difference in the results.
   txt = ['_D' int2str(10*Dia)                              ...
          '_V'  int2str(100*Vinf) '_I' int2str(1000*I)      ...
          '_Lu' int2str(Lu)       '_df' int2str(10000*df)   ...
          '_nf' int2str(nf)       '_W'  int2str(Wop*1000) '.bin'];
   save('-float-binary',['Spsi' txt],'Spsi');

   clear Spsi;

   if (flagg ~= 0)
      % The derivative wrt the selected parameter.  Use duplicate 
      % variables to avoid memory limitations.  Note that SpsiI
      % has been redefined such that it is the real and imaginary
      % part of the complex step perturbation.
      dS    = zeros(size(SpsiI),'single');
      dS(:) = SpsiI/epsilon;
      save('-float-binary',['dSd' ch txt],'dS');
      clear dS SpsiI;
   end

   if (towerFlag == 1)

      % Find the transforms which may be dependent on Vinf.
      [Ty_g,Td_n,Tp_b] = buildJointTransforms (chi0,psi0,beta0);
      Tg_r = (Tn_y.')*(Ty_g.');

      sgr_g = Ogy_g + Ty_g*(Oyd_y + Tn_y*Td_n*[0;0;Odp_d(3)]);

      [SpsiTR,SpsiTI] = mbcSijTower (Neb,Net,nf,df,Vinf,I,Lu,Wnom, ...
                                     r3,zet_g,sgr_g,Tg_r,Ty_g,psi0);

      SpsiTower = real(SpsiTR) + i*real(SpsiTI);
      save('-float-binary',['SpsiTower' txt],'SpsiTower');
      clear SpsiTower;

      if (flagg ~= 0)
         dSTower = (imag(SpsiTR) + i*imag(SpsiTI))/epsilon;
         save('-float-binary',['dSTowerd' ch txt],'dSTower');
         clear dSTower;
      end

   end

end


