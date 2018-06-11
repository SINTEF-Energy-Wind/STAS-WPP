function [eta,] = aeroelasticOperatingPoint ...
                     ()
%
% Compute a steady-state operating point using Newton's method.
%
% Version:        Changes:
% --------        -------------
% 06.06.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 06.06.2018      
%
% Inputs:
% -------
% cnv             : Convergence criterion.
% Ns              : Maximum number of Newton iterations.
% bfac            : Step factor for each iteration.  It is wise to
%                   force a couple small steps at the beginning, so
%                   that the initial step doesn't put the estimate
%                   at a point far from the correct solution, where
%                   nonlinear features may prevent convergence.
%                   Something like 0.1,0.2,0.4,0.5,1.0,1.0... has
%                   been observed to work.
%
% Outputs:
% --------
%                 : 

% Lock the foundation DOFs, if there are no soil springs; and lock
% all the joint DOFs at their prescribed values.
ldofs = [fndlock;Neta-6+[1:6].'];
deldofs = [ldofs;Neta+ldofs];

% Do an initial calculation at the guess state to get an initial
% Rval.
[Lmat,Rvec,y] = aeroelasticNL (0,etasa0,y1,Fext,grav,    ...
                               psiFlag,s,a,P,ret,slv,    ...
                               shape,mdamp,              ...
                               Tas,Try,Vg,ch,Lel,foilwt, ...
                               aoaz,aoast,xas,yas,Psi);
Rvpsi = TBpsi_Nst*Rvec ...
      - Waero*TBpsi_Nst*Lmat*dTpsiB_Nst*TBpsi_Nst*etasa0;
[Rvpar,dret,jnkc] = partitionMatrix(Rvpsi,deldofs,[]);
Ndr  = size(dret,1);
Rvr  = Rvpar(1:Ndr);
Rval = (Rvr.')*Rvr;

q   = y(1:Ndj);      % These now agree with the initial x.
dq  = y(Ndj+[1:Ndj]);
d2q = zeros(Ndj,1);  % Only used for structureLin call.  Want this
                     % Want this to be zero to match Rvec
F      = y(3*Ndj+6*Nel+[1:Ndj]);
etapsi = TBpsi_Nst*etasa0;

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
   F   = y(3*Ndj+6*Nel+[1:Ndj]);

   % ---------------------------------------------------------------
   % The first stage is to get an update equation for x.
   [LHS,A,Bu,By,C,Du,Dy,shape,freq,mdamp,jnkr,jnks,Psi] =  ...
               aeroelasticLin (eta,Vg,s,a,q,dq,d2q,P,F,    ...
                               psiFlag,modeflag,shpflag,   ...
                               shape0,mdamp0,grav);
%return
   % Link the linear model.

   % Large D matrix, super slow!
%   Dnorm = ones(Nylin,1);
%   [AA,BB,CC,DD] = modularToUnifiedStateSpace (A,Bu,By,C,Du,Dy,Dnorm);

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
      [Lmat,Rvec,y] = aeroelasticNL (t0,etaB,y1,Fext,grav,     ...
                                     psiFlag,s,a,P,ret,slv,    ...
                                     shape0,mdamp0,            ...
                                     Tas,Try,Vg,ch,Lel,foilwt, ...
                                     aoaz,aoast,xas,yas,Psi);

      % Check acceptability of Rvec.
      Rvpsi = TBpsi_Nst*Rvec - Waero*TBpsi_Nst*Lmat*dTpsiB_Nst*etapsi1;
      [Rvpar,dret,jnkc] = partitionMatrix(Rvpsi,deldofs,[]);
      Rvr = Rvpar(1:Ndr);
      R1 = (Rvr.')*Rvr;

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
F   = y(3*Ndj+6*Nel+[1:Ndj]);
Tgen = y(size(y,1));

[LHS,A,Bu,By,C,Du,Dy,shape,freq,mdamp,jnkr,jnks,Psi] =  ...
            aeroelasticLin (eta,Vg,s,a,q,dq,d2q,P,F,    ...
                            psiFlag,modeflag,shpflag,   ...
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


