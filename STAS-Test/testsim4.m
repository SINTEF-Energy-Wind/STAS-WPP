% Test program for a Newton's method solution of the steady-state
% operating point.

clear;

rand('state',100);

%[s,a] = STASTurbine_aeroValidation ();
%[s,a] = STASTurbine_DTU10MW ();
[s,a] = STASTurbine_NREL5MW_land ();
%[s,a] = STASTurbine_Tjaereborg ();

a.icp = 0;  % Force full aero states.

[length,time,mass,current,velocity,force,power,stress, ...
 ndens,nvisc,stiffness,damping,resistance,inductance,  ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
[imdofs,Nmd] = getmdofRefs (s);

Nb = s.Nb;
Neb = s.blade(1).Nel;
Nel = Nb*Neb;

% For linear analysis, the azimuth angle is arbitrary, as its
% influence will be eliminated by the MBC transform.
azi0 = 90*pi/180;
cp0 = cos(azi0);
sp0 = sin(azi0);

% ===================================================================
% Input parameters defining the load case.
a.dens = 1.225/(mass/(length^3));
a.visc = 1.789e-5/(force*time/(length^2));

betas = 0.0*ones(3,1)*pi/180;
Vmag = 8/velocity;
yaw = 0*pi/180;
Vinf = [Vmag*cos(yaw);Vmag*sin(yaw);0];
Omega = 0.9524*time; % 12.0*(2*pi/60)*time;
alf = 0.17;
Dt = 4.5*ones(Nel,1)/length;
grav = [0;0;-9.807]/(length/(time^2));

t0 = 0;

% ===================================================================
% Best if psiFlag = 1, as this gives the appropriate induced velocity
% in yawed flow conditions.
psiFlag  = 1;
modeflag = 1;
shpflag  = 0;

% Fill out the global velocity vector based on the inputs.
Vg = zeros(3*Nel,1);
Vg(1:3:3*Nel-2) = Vinf(1);
Vg(2:3:3*Nel-1) = Vinf(2);
Vg(3:3:3*Nel)   = Vinf(3);

% Initial guess for induced velocity.
Vi0 = [-0.35*Vinf(1);0];
Viguess = zeros(2*Nel,1);
Viguess(1:2:2*Nel-1) = Vi0(1);
Viguess(2:2:2*Nel)   = Vi0(2);

% ===================================================================
% Initialize the structral calculation.
[xs0,etas0,q0,dq0dt,d2q0dt2,P,shape0,shs,freq0,mdamp0,ret,slv, ...
 Ndj,Nnod,Neta] = structInit (s,yaw,azi0,betas,Omega,modeflag);

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);
Nret = size(ret,1);
Nslv = size(slv,1);
Nylin = 4*Ndj + 9*Nnod + 137*Nel + 5 + 1;  % + 1: Tgen
Nu = Ndj + 3*Nel;

Ydof  = idofs(3);
Ddof  = idofs(4);
nodof = idofm(6) - 6;

iazi = Nylin - 5;
iW   = Nylin - 4;

% ===================================================================
% Eliminate DOFs from the nonlinear and linear models.  For the
% present study, lock the joints except for the driveshaft, and get
% rid of the azimuth, which isn't needed with MBC.  Also lock the
% foundation base DOFs, if this is needed.
%deldofs = [Neta-6+[1 2 3 4 5 6] 2*Neta-6+[1 3 4 5 6]].';
%deldofs = [Neta-6+[1 3 4 5 6] 2*Neta-6+[1 3 4 5 6]].';
deldofs = [[1:6] Neta-6+[1 2 3 4 5 6] Neta+[1:6] 2*Neta-6+[1 3 4 5 6]].';

% ===================================================================
% BEM setup, parameters that stay fixed throughout the calculations,
% and initial estimates of the aerodynamic states.
[Tas,ch,Lel,foilwt,aoaz,aoast,xas,yas,iq] = BEMsetup (s,a);

Td_n = [cp0 -sp0 0;sp0 cp0 0;0 0 1];
Try = Tn_y;
Tyy0 = TFromTheta (q0(idofs(3)+[4:6]));
Ty0g = TFromTheta (P(idofs(3)+[4:6]));
Trg = Ty0g*Tyy0*Try;

[Tar,Trg,TB0g,TsB,TBB0,dTar,dTsB,dTBB0,wg] = ...
                      BEMprepTransforms (s,a,q0,dq0dt,P,Tas);

[zr,Area,Dp,rp,Lp,xeg,xhg,xyg] = ...
         BEMprepProjections (s,a,q0,dq0dt,P,Try,Trg);
Dps = zeros(Nel,1);
Dps(1:Neb)         = Dp(1);
Dps(Neb+[1:Neb])   = Dp(2);
Dps(2*Neb+[1:Neb]) = Dp(3);

[Waero,azi,Dy] = rotorSpeedAero (q0,dq0dt,P,Try,Ydof,Ddof,nodof);
xa0 = BEMinit (psiFlag,                         ...
               Viguess,Tar,Trg,Vg,wg,zr,ch,Lel, ...
               aoast,a.dens,Area,Dps,azi,Waero);

% Best estimate of the initial modal aero states.
bsh = bladeModeShape (s,ret,shape0);
Psi = aeroPsi (a,rp,bsh);
etaa0 = (Psi.')*xa0;

% ===================================================================
% Prepare MBC transforms.  These are needed for the reduced set of 
% structural and aero states.  (Presently there is no pitch control,
% so I don't need to transform the controls.)
if (a.icp(1) == 0)
   Ncp  = Neb;
   Nsta = size(etaa0,1);
else
   Ncp  = size(a.icp,1);
   Nsta = 7*Ncp*a.Nb;
end
Nst = 2*Neta + Nsta;
if (modeflag == 0)
   [b1_Ns,b2_Ns,b3_Ns]       = MBCindices_Nret (Ndj,idofs,ret,slv);
elseif (modeflag == 1)
   [b1_Ns,b2_Ns,b3_Ns]       = MBCindices_Nmd (Neta,imdofs);
end
[b1_Nxa,b2_Nxa,b3_Nxa]       = MBCindices_N (7*Ncp);

b1_Nst = [b1_Ns;Neta+b1_Ns;2*Neta+b1_Nxa];
b2_Nst = [b2_Ns;Neta+b2_Ns;2*Neta+b2_Nxa];
b3_Nst = [b3_Ns;Neta+b3_Ns;2*Neta+b3_Nxa];
[TpsiB_Nst,TBpsi_Nst]     = MBC (Nst,b1_Nst,b2_Nst,b3_Nst,azi0);
[dTpsiB_Nst,dTBpsi_Nst]   = derivMBC (Nst,b1_Nst,b2_Nst,b3_Nst,azi0);
[d2TpsiB_Nst,d2TBpsi_Nst] = secondDerivMBC (Nst,b1_Nst,b2_Nst,b3_Nst,azi0);

% Uncomment this to deactivate MBC.
%TpsiB_Nst = speye(Nst);
%TBpsi_Nst = speye(Nst);
%dTpsiB_Nst = sparse(Nst,Nst);
%dTBpsi_Nst = sparse(Nst,Nst);
%d2TpsiB_Nst = sparse(Nst,Nst);
%d2TBpsi_Nst = sparse(Nst,Nst);


% ===================================================================
% Get initial values for Rval, q, dq/dt.
etasa0 = [etas0;etaa0];
[q,dq] = qFromx (etas0,s,P,shape0,ret,slv,q0(slv));
y1 = [q;dq;zeros(2*Ndj+6*Nel+1,1)];
[Tar,Trg,TB0g,TsB,TBB0,dTar,dTsB,dTBB0,wg] = ...
              BEMprepTransforms (s,a,q,dq,P,Tas);
[zr,Area,Dp,rp,Lp,xeg,xhg,xyg] = ...
              BEMprepProjections (s,a,q,dq,P,Try,Trg);

[Lmat,Rvec,y] = turbineNL_trial (t0,etasa0,y1,             ...
                                 psiFlag,s,a,P,ret,slv,    ...
                                 shape0,mdamp0,grav,       ...
                                 Tas,Try,Vg,ch,Lel,foilwt, ...
                                 aoaz,aoast,xas,yas,Psi);
Rvpsi = TBpsi_Nst*Rvec ...
      - Waero*TBpsi_Nst*Lmat*dTpsiB_Nst*TBpsi_Nst*etasa0;
[Rvpar,dret,jnkc] = partitionMatrix(Rvpsi,deldofs,[]);
Ndr  = size(dret,1);
Rvr  = Rvpar(1:Ndr);
Rval = (Rvr.')*Rvr;

Rval

q   = y(1:Ndj);      % These now agree with the initial x.
dq  = y(Ndj+[1:Ndj]);
d2q = zeros(Ndj,1);  % y(2*Ndj+[1:Ndj]); % Only used for structureLin call.
                                         % Want this to be zero to match
                                         % Rvec.
F      = y(3*Ndj+6*Nel+[1:Ndj]);
etapsi = TBpsi_Nst*etasa0;

% ===================================================================
% Newton's method solution for the steady state.

tic;

meth  = 1;   % 1: linearized eqs.  2: Broyden.
cnv   = eps^(0.6);
Ns    = 10; % 100;
bta   = [0.1 0.2 0.4 0.5 0.5 1.0*ones(1,97)].';
conv  = 0;
iter  = 0;
Niter = 0;
while ((Rval > cnv) && (iter < Ns))
   iter = iter + 1;

'================================='
iter

   eta = TpsiB_Nst*etapsi;

   q   = y(1:Ndj);
   dq  = y(Ndj+[1:Ndj]);
   d2q = zeros(Ndj,1);
%d2q = y(2*Ndj+[1:Ndj]);
   F   = y(3*Ndj+6*Nel+[1:Ndj]);

   % ---------------------------------------------------------------
   % The first stage is to get an update equation for x.  This can
   % use either the tangent dynamics or Broyden (secant).
   if (meth == 1)

      [LHS,A,Bu,By,C,Du,Dy,shape,freq,mdamp,jnkr,jnks,Psi] =  ...
                  turbineLin_trial (eta,Vg,s,a,q,dq,d2q,P,F,  ...
                                    psiFlag,modeflag,shpflag, ...
                                    shape0,mdamp0,grav);
%return
      % Link the linear model.

      % Large D matrix, super slow!
%      Dnorm = ones(Nylin,1);
%      [AA,BB,CC,DD] = modularToUnifiedStateSpace (A,Bu,By,C,Du,Dy,Dnorm);

      % Use a procedure that takes advantage of the structure
      % of Dy... super fast!
      adofs = 4*Ndj+9*Nnod+[1:137*Nel].';
      rdofs = [[1:4*Ndj+9*Nnod] [Nylin-5:Nylin]].';  % Incl. Tgen.
      Na = size(adofs,1);
      Nr = size(rdofs,1);
      Nz = size(C,2) + size(Du,2);
      H = speye(Nylin) - Dy;
      Hrr = H(rdofs,rdofs);
      Hra = H(rdofs,adofs);
      Har = H(adofs,rdofs);
      Haa = H(adofs,adofs);
      R   = [C Du];
      Ra  = R(adofs,:);
      Rr  = R(rdofs,:);
      HHR = Haa\[Har Ra];
      H4  = Hrr - Hra*HHR(:,1:Nr);
      Sr  = H4\(Rr - Hra*HHR(:,Nr+[1:Nz]));
      Sa  = -HHR(:,1:Nr)*Sr + HHR(:,Nr+[1:Nz]);
      CC  = sparse(Nylin,Nst);
      DD  = sparse(Nylin,Nu);
      CC(rdofs,:) = Sr(:,1:Nst);
      CC(adofs,:) = Sa(:,1:Nst);
      DD(rdofs,:) = Sr(:,Nst+[1:Nu]);
      DD(adofs,:) = Sa(:,Nst+[1:Nu]);
      AA = A + By*CC;
      BB = Bu + By*DD;

%   else
      % Broyden.
      
   end

   % MBC transform and DOF elimination of the A matrix.
   [Waero,azi,jnk] = rotorSpeedAero (q,dq,P,Try,Ydof,Ddof,nodof);
   Apsi = TBpsi_Nst*AA*TpsiB_Nst - Waero*TBpsi_Nst*LHS*dTpsiB_Nst;
   Bypsi = spalloc (Nst,Nylin,2*Nylin);
   Bypsi(:,iW) = -TBpsi_Nst*LHS*dTpsiB_Nst*etapsi;
   AApsi = Apsi + Bypsi*CC*TpsiB_Nst;
   [Appar,jnkr,jnkc] = partitionMatrix(AApsi,deldofs,deldofs);
   Ar  = Appar(1:Ndr,1:Ndr);

   % ----------------------------------------------------------------
   % Apply Newton's method on the MBC-transformed equations.
   lam = 1;
   deta = -Ar\Rvr;
   lflg = 0;
   litmax = 10;
   liter = 0;
   while ((lflg == 0) && (liter < litmax))
      liter = liter + 1;

      % Try a step.  Include e.g. constant shaft speed, that is not
      % influenced by the dret DOFs.
      etapsi1 = TBpsi_Nst*etasa0; 
      etapsi1(dret) = etapsi(dret) + bta(iter)*lam*deta;
      etaB = TpsiB_Nst*etapsi1;

      % Evaluate the step.
      %
      % We now need to update y to match x.  The first thing is that the
      % q's must agree with the x's, as the structural positions and
      % velocities are needed in the computation of aerodynamic forces.
      [q,dq] = qFromx (etaB(1:2*Neta),s,P,shape0,ret,slv,q0(slv));
      [Waero,azi,jnk] = rotorSpeedAero (q,dq,P,Try,Ydof,Ddof,nodof);

      % Accelerations are not used in the call to turbineNL.  Likewise,
      % airfoil forces will be updated from within turbineNL.  These 
      % can for the time being be zero.
      y1 = [q;dq;zeros(2*Ndj+6*Nel+1,1)];

      % Aero basis functions.
      [Tar,Trg,TB0g,TsB,TBB0,dTar,dTsB,dTBB0,wg] = ...
                    BEMprepTransforms (s,a,q,dq,P,Tas);
      [zr,Area,Dp,rp,Lp,xeg,xhg,xyg] = ...
                    BEMprepProjections (s,a,q,dq,P,Try,Trg);
      bsh = bladeModeShape (s,ret,shape0);
      Psi = aeroPsi (a,rp,bsh);

      % Nonlinear equations.
      [Lmat,Rvec,y] = turbineNL_trial (t0,etaB,y1,               ...
                                       psiFlag,s,a,P,ret,slv,    ...
                                       shape0,mdamp0,grav,       ...
                                       Tas,Try,Vg,ch,Lel,foilwt, ...
                                       aoaz,aoast,xas,yas,Psi);

      % Check acceptability of Rvec.
      Rvpsi = TBpsi_Nst*Rvec - Waero*TBpsi_Nst*Lmat*dTpsiB_Nst*etapsi1;
      [Rvpar,dret,jnkc] = partitionMatrix(Rvpsi,deldofs,[]);
      Rvr = Rvpar(1:Ndr);
      R1 = (Rvr.')*Rvr;



%[[[1:375] [1:375] [1:252]].' Rvec]



full([liter Rval R1])

val = toc();
fid = fopen('stat.txt','a');
fprintf(fid,'%6d %10.4f %10.4f %9.4f %9.4f %9.4f %7.1f\n', ...
        iter,R1,Rval,F(idofs(6)+8*6+3),etaB(imdofs(6)+1),  ...
        etaB(2*Neta+1),val);
fclose(fid);

      if (R1 < Rval)
         % OK!  Prepare for the next iteration.
         lflg = 1;
         Rval = R1;
         etapsi = etapsi1;

      else
         % Backtrack.
         lam = 0.5*lam;
         if (liter == litmax)
            iter
            'Warning, proceeding without lambda convergence.'
            lflg = 1;
            Rval = R1;
            etapsi = etapsi1;

         end
      end

   end

[toc Rval]

val = toc();
fid = fopen('stat.txt','a');
fprintf(fid,'%6d %+5.6e %10.6f %10.6f %10.3f\n', ...
        iter,Rval,F(idofs(6)+8*6+3),etaB(imdofs(6)+1),val);
fprintf(fid,'--------------------------\n');
fclose(fid);

end

% Final matrices.
eta = TpsiB_Nst*etapsi;

q   = y(1:Ndj);
dq  = y(Ndj+[1:Ndj]);
d2q = zeros(Ndj,1);
%d2q = y(2*Ndj+[1:Ndj]);
F   = y(3*Ndj+6*Nel+[1:Ndj]);
Tgen = y(size(y,1));

[LHS,A,Bu,By,C,Du,Dy,shape,freq,mdamp,jnkr,jnks,Psi] =  ...
            turbineLin_trial (eta,Vg,s,a,q,dq,d2q,P,F,  ...
                              psiFlag,modeflag,shpflag, ...
                              shape0,mdamp0,grav);
adofs = 4*Ndj+9*Nnod+[1:137*Nel].';
rdofs = [[1:4*Ndj+9*Nnod] [Nylin-5:Nylin]].';
Na = size(adofs,1);
Nr = size(rdofs,1);
Nz = size(C,2) + size(Du,2);
H = speye(Nylin) - Dy;
Hrr = H(rdofs,rdofs);
Hra = H(rdofs,adofs);
Har = H(adofs,rdofs);
Haa = H(adofs,adofs);
R   = [C Du];
Ra  = R(adofs,:);
Rr  = R(rdofs,:);
HHR = Haa\[Har Ra];
H4  = Hrr - Hra*HHR(:,1:Nr);
Sr  = H4\(Rr - Hra*HHR(:,Nr+[1:Nz]));
Sa  = -HHR(:,1:Nr)*Sr + HHR(:,Nr+[1:Nz]);
CC  = sparse(Nylin,Nst);
DD  = sparse(Nylin,Nu);
CC(rdofs,:) = Sr(:,1:Nst);
CC(adofs,:) = Sa(:,1:Nst);
DD(rdofs,:) = Sr(:,Nst+[1:Nu]);
DD(adofs,:) = Sa(:,Nst+[1:Nu]);
AA = A + By*CC;
BB = Bu + By*DD;

[Waero,azi,jnk] = rotorSpeedAero (q,dq,P,Try,Ydof,Ddof,nodof);
Apsi = TBpsi_Nst*AA*TpsiB_Nst - Waero*TBpsi_Nst*LHS*dTpsiB_Nst;
Bypsi = spalloc (Nst,Nylin,2*Nylin);
Bypsi(:,iW)   = -TBpsi_Nst*LHS*dTpsiB_Nst*etapsi;
AApsi = Apsi + Bypsi*CC*TpsiB_Nst;
Lpsi = TBpsi_Nst*LHS*TpsiB_Nst;

[Appar,jnkr,jnkc] = partitionMatrix(AApsi,deldofs,deldofs);
[Lpar,jnkr,jnkc]  = partitionMatrix(Lpsi,deldofs,deldofs);
Ar  = Appar(1:Ndr,1:Ndr);
Lr  = Lpar(1:Ndr,1:Ndr);

[slap,shp,ifrq] = eigVal (Lr\Ar);



% ========================================================
% Model for Guodong.

% Sum over Vg's to get a single rotor-avg V input to the model.
BV = sum(BB(:,Ndj+[1:3:3*Nel]),2);
BVpsi = TBpsi_Nst*BV;

% Lock only the joint DOFs, let the foundation ref node stay.
ddofs = [Neta-6+[1 2 3 4 5 6] 2*Neta-6+[1 3 4 5 6]].';

[Appar,dret2,jnkc] = partitionMatrix(AApsi,ddofs,ddofs);
[Lpar,jnkr,jnkc]   = partitionMatrix(Lpsi,ddofs,ddofs);
[Bpar,jnkr,jnkc]   = partitionMatrix(BVpsi,ddofs,[]);
Ndr2 = size(dret2,1);
Ar  = Appar(1:Ndr2,1:Ndr2);
Lr  = Lpar(1:Ndr2,1:Ndr2);
Br  = Bpar(1:Ndr2,:);

qn1 = zeros(6,1);
qn2 = q(7:12);
Pn2 = P(7:12);
Pn1 = zeros(6,1);
Pn1(4:6) = Pn2(4:6);
[xe,TsB] = elementCSFromNodes (qn1,qn2,Pn1,Pn2);
dTsB = derivElementCS (qn1,qn2,Pn1,Pn2,TsB);
dmu = dmudq (qn1,qn2,Pn1,Pn2,TsB,dTsB);
kes = buildkes (s.foundation.EEs(:,1:6),s.foundation.Lel(1));
Cr = kes(6,:)*dmu*[sparse(6,Ndr2); ...
                   [shape(7:12,1:Neta-6) sparse(6,Neta-6+1+(7*48))]];

fidA = fopen('A.txt','w');
fidN = fopen('N.txt','w');
fidB = fopen('B.txt','w');
fidC = fopen('C.txt','w');
fprintf(fidA,'%10d %10d\n',size(Ar,1),size(Ar,2));
fprintf(fidN,'%10d %10d\n',size(Lr,1),size(Lr,2));
fprintf(fidB,'%10d %10d\n',size(Br,1),size(Br,2));
fprintf(fidC,'%10d %10d\n',size(Cr,1),size(Cr,2));
for icol = 1:Ndr2
   for irow = 1:Ndr2
      fprintf(fidA,'%+5.6e\n',Ar(irow,icol));
      fprintf(fidN,'%+5.6e\n',Lr(irow,icol));
   end
   fprintf(fidB,'%+5.6e\n',Br(icol));
   fprintf(fidC,'%+5.6e\n',Cr(icol));
end

fclose('all');

%save sim4vars.txt
