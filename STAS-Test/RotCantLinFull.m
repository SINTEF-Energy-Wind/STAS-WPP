function [M,C,K] = RotCantLinFull (q,dqdt,d2qdt2,P,mes,kes,adamp)
%
% q consists of qB, then qn's.
%

Ndof = size(q,1);
Nel = (Ndof - 6)/6;

r = sqrt(P(1)^2 + P(2)^2);

Mf = zeros(Ndof,Ndof);  % Full matrices including qB DOFs.
Cf = zeros(Ndof,Ndof);
Kf = zeros(Ndof,Ndof);
for iel = 1:Nel

   ic6 = 6*(iel-1);

   if (iel == 1)
      qn1      = zeros(6,1);
      dqn1dt   = zeros(6,1);
      d2qn1dt2 = zeros(6,1);
      Pn1      = zeros(6,1);
      Pn1(4:6) = P(10:12);
   else
      qn1      = q(ic6+[1:6]);
      dqn1dt   = dqdt(ic6+[1:6]);
      d2qn1dt2 = d2qdt2(ic6+[1:6]);
      Pn1      = P(ic6+[1:6]);
   end

   qB       = q(1:6);
   dqBdt    = dqdt(1:6);
   d2qBdt2  = d2qdt2(1:6);
   PB       = P(1:6);
   qn2      = q(ic6+[7:12]);
   dqn2dt   = dqdt(ic6+[7:12]);
   d2qn2dt2 = d2qdt2(ic6+[7:12]);
   Pn2      = P(ic6+[7:12]);

   dqedt   = [dqBdt;dqn1dt;dqn2dt];
   d2qedt2 = [d2qBdt2;d2qn1dt2;d2qn2dt2];

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

   d2qG = zeros(18,1);
   d2qL = d2qedt2;

   [mme,dmmeL,dmmeG,gge,dgge,dgged,hv,dhhe,dhhed,kv,dkke] = ...
              buildElementLin (mes,kes,dqedt,d2qL,d2qG,       ...
                               Qu1,Qu2,dQu1,dQu2,d2Qu1,d2Qu2, ...
                               mu,dmu,d2mu,TsB,dTsB,d2TsB);

   [cde,kde] = dampingCLin (adamp*kes,dqedt,dmu,d2mu);

   Mf(1:6,1:6) = Mf(1:6,1:6) + mme(1:6,1:6);
   Cf(1:6,1:6) = Cf(1:6,1:6) + gge(1:6,1:6) + dgged(1:6,1:6) - dhhed(1:6,1:6) + cde(1:6,1:6);
   Kf(1:6,1:6) = Kf(1:6,1:6) - dmmeL(1:6,1:6) - dmmeG(1:6,1:6) + dgge(1:6,1:6) ...
                             - dhhe(1:6,1:6)  + dkke(1:6,1:6)  + kde(1:6,1:6);

   if (iel == 1)
      % Elastic DOFs 7:12 in mme,cce,kke do not exist in Mf,Cf,Kf.
      Mf(1:6,ic6+[7:12]) = Mf(1:6,ic6+[7:12]) + mme(1:6,13:18);
      Mf(ic6+[7:12],1:6) = Mf(ic6+[7:12],1:6) + mme(13:18,1:6);
      Mf(ic6+[7:12],ic6+[7:12]) = Mf(ic6+[7:12],ic6+[7:12]) ...
                                + mme(13:18,13:18);
      Cf(1:6,ic6+[7:12]) = Cf(1:6,ic6+[7:12])                ...
                         + gge(1:6,13:18) + dgged(1:6,13:18) ...
                         - dhhed(1:6,13:18) + cde(1:6,13:18);
      Cf(ic6+[7:12],1:6) = Cf(ic6+[7:12],1:6)                ...
                         + gge(13:18,1:6) + dgged(13:18,1:6) ...
                         - dhhed(13:18,1:6) + cde(13:18,1:6);
      Cf(ic6+[7:12],ic6+[7:12]) = Cf(ic6+[7:12],ic6+[7:12])      ...
                         + gge(13:18,13:18) + dgged(13:18,13:18) ...
                         - dhhed(13:18,13:18) + cde(13:18,13:18);
      Kf(1:6,ic6+[7:12]) = Kf(1:6,ic6+[7:12])                                    ...
                         - dmmeL(1:6,13:18) - dmmeG(1:6,13:18) + dgge(1:6,13:18) ...
                         - dhhe(1:6,13:18)  + dkke(1:6,13:18)  + kde(1:6,13:18);
      Kf(ic6+[7:12],1:6) = Kf(ic6+[7:12],1:6)                                    ...
                         - dmmeL(13:18,1:6) - dmmeG(13:18,1:6) + dgge(13:18,1:6) ...
                         - dhhe(13:18,1:6)  + dkke(13:18,1:6)  + kde(13:18,1:6);
      Kf(ic6+[7:12],ic6+[7:12]) = Kf(ic6+[7:12],ic6+[7:12])                            ...
                         - dmmeL(13:18,13:18) - dmmeG(13:18,13:18) + dgge(13:18,13:18) ...
                         - dhhe(13:18,13:18)  + dkke(13:18,13:18)  + kde(13:18,13:18);

   else
      Mf(1:6,ic6+[1:12]) = Mf(1:6,ic6+[1:12]) + mme(1:6,7:18);
      Mf(ic6+[1:12],1:6) = Mf(ic6+[1:12],1:6) + mme(7:18,1:6);
      Mf(ic6+[1:12],ic6+[1:12]) = Mf(ic6+[1:12],ic6+[1:12]) ...
                                + mme(7:18,7:18);
      Cf(1:6,ic6+[1:12]) = Cf(1:6,ic6+[1:12])              ...
                         + gge(1:6,7:18) + dgged(1:6,7:18) ...
                         - dhhed(1:6,7:18) + cde(1:6,7:18);
      Cf(ic6+[1:12],1:6) = Cf(ic6+[1:12],1:6)              ...
                         + gge(7:18,1:6) + dgged(7:18,1:6) ...
                         - dhhed(7:18,1:6) + cde(7:18,1:6);
      Cf(ic6+[1:12],ic6+[1:12]) = Cf(ic6+[1:12],ic6+[1:12])  ...
                         + gge(7:18,7:18) + dgged(7:18,7:18) ...
                         - dhhed(7:18,7:18) + cde(7:18,7:18);
      Kf(1:6,ic6+[1:12]) = Kf(1:6,ic6+[1:12])                                 ...
                         - dmmeL(1:6,7:18) - dmmeG(1:6,7:18) + dgge(1:6,7:18) ...
                         - dhhe(1:6,7:18)  + dkke(1:6,7:18)  + kde(1:6,7:18);
      Kf(ic6+[1:12],1:6) = Kf(ic6+[1:12],1:6)                                 ...
                         - dmmeL(7:18,1:6) - dmmeG(7:18,1:6) + dgge(7:18,1:6) ...
                         - dhhe(7:18,1:6)  + dkke(7:18,1:6)  + kde(7:18,1:6);
      Kf(ic6+[1:12],ic6+[1:12]) = Kf(ic6+[1:12],ic6+[1:12])                      ...
                         - dmmeL(7:18,7:18) - dmmeG(7:18,7:18) + dgge(7:18,7:18) ...
                         - dhhe(7:18,7:18)  + dkke(7:18,7:18)  + kde(7:18,7:18);
   end

end

% Constrain the root DOFs.
M = Mf(7:Ndof,7:Ndof);
C = Cf(7:Ndof,7:Ndof);
K = Kf(7:Ndof,7:Ndof);

