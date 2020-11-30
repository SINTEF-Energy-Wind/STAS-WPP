function [Lpsi,Rpsi,ypsi,Apsi,Bpsi,Cpsi,Dpsi] =                   ...
                  MBCCLT (linFlag,xpsi,dxpsi,upsi,s,a,            ...
                          epar,ppar,ypar,c,grav,P,shape0,mdamp0,  ...
                          Tas,Try,ch,Lel,foilwt,aoaz,aoast,       ...
                          xas,yas,Psi0,igen,ipit,iyaw)
%
% Return the multi-blade coordinate transformed version of the dynamic
% closed-loop turbine equations.
%
% Note that d2eta/dt2 is an element of the input dxpsi, and is therefore
% excluded from the upsi vector; however, the indexing of the upsi
% DOFs in the matrices include the d2eta/dt2 entries.
%
%   States:              y vector:                u vector:
%                                                 F       Ndj   (Env)
%                                              d2eta/dt2  Neta  (sol'n)
%                                                 Vg    3*Nae   (Env)
%                                                 we        1   (Grid)
%                                                 th_e      2   (Grid)
%                                                 vsd,q   3,4   (Grid)
%                                                 Pc        5   (Pl.cont.)
%                                                 Qc        6   (Pl.cont.)
% ----------------------- Structure --------------------------
%   eta        N         q         Ndj            F       Ndj   (u)
%   deta/dt    N         dq/dt     Ndj         d2eta/dt2  Neta  (u) (NOT in upsi)
%                        d2q/dt2   Ndj
%                        F         Ndj
%
% ---------------------- Aerodynamic -------------------------
%   ad         1                                  Vg        3   (u)
%   a1,a2     2:3
%   Vih z,t   4:5
%   Vi z,t    6:7
%
% ------------------------ Actuators --------------------------
% (Repeat for each blade and yaw system.)
%   ba         1         b,dbdt    (-)  (Turb)     bhat      1   (Cont)
%   dbadt      2         Ta         1
%
% ----------------------- Electrical --------------------------
%   igd,q     1,2        wg        (-)  (Turb)     we        1   (u)
%   Vdc        3         Tg         1              th_e      2   (u)
%   ipd,q     4,5        isd,q     2,3             vsd,q    3,4  (u)
%   imgd,q    6,7                                  ihgd,q   5,6  (Cont)
%   Psig      8,9                                  Vhdc      7   (cpar)
%   wemg      10                                   Qh        8   (u)
%   th_m      11
%   vmsd,q   12,13
%   Psie      14
%   impd,q   15,16
%   imsd,q   17,18
%   vmsd,q   19,20
%   Vmdc      21
%   PsiDC     22
%   PsiQ      23
%   PsiP     24,25
%
% ------------------------- Control ---------------------------
%   W*         1         bhat      1:3             Pc        1   (u)
%   V*         2         ihgd,q    4:5             azi       2   (turbine)
%   Wm         3         yhat       6              W         3   (turbine)
%   betm       4                                   bet      4:6  (turbine)
%   PsiWb      5                                   Pe        7   (turbine)
%   Pf         6                                   vT       8:9  (turbine)
%   PehRSC     7                                   Mbl     10:12 (turbine)(moment or th_nod)
%   ntw       8,9                                  wang     13   (turbine)(derived?)
%   ntb      10,11
%   blp       12
%   wgm       13  
%   iwg       14
%   Pem       15
%   PsiP      16
%   nd       17,18
%   vmF       19
%   ivmF      20
%   vdF       21
%   vmS       22
%   ivmS      23
%   vdS       24
%   vmS       25
%   ivmS      26
%   vdS       27
%   Mpsim    28,29
%   PsiM     30,31
%
% Version:        Changes:
% --------        -------------
% 14.02.2019      Original code.
% 23.01.2020      Updated to accommodate revised equations-of-motion
%                 and u vector.
%
% Version:        Verification:
% --------        -------------
% 14.02.2019      
% 23.01.2020      A,C matrix derivatives verified by complex step
%                 with L,R outputs.
%

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
Ndj  = Ndof + 6;
[imdofs,Neta] = getmdofRefs (s);
Neb  = s.blade(1).Nel;
Nxa  = size(Psi0,2);
Nxp  = 6;
Nx1  = 2*Neta + Nxa;
Nx   = size(xpsi,1);

Nun  = size(upsi,1);
Nu   = Nun + Neta;
Ny   = 4*Ndj + 4 + 3 + 6;

iW   = 2*Neta - 4;
iazi = Neta - 4;
azi  = xpsi(iazi);

[blxdof,bludof,blydof] = MBCindices (Nxa,Nxp,Neta,Ndj,Neb,imdofs,idofs);

[x,Txpx,dTxpx]  = buildTxpx (xpsi,dxpsi,blxdof(:,1),blxdof(:,2),blxdof(:,3),iazi);
[TpBx,TBpx]     = MBC (Nx,blxdof(:,1),blxdof(:,2),blxdof(:,3),azi);
[TpBu,TBpu]     = MBC (Nu,bludof(:,1),bludof(:,2),bludof(:,3),azi);
[TpBy,TBpy]     = MBC (Ny,blydof(:,1),blydof(:,2),blydof(:,3),azi);

% Don't neet extra azimuth-column terms or complex-step gradients
% in the premultiplying transform.  This stems from the fact that
% the azimuth and rotor speed states are the same in body and MBC
% coordinates.  So the state equation dPsi/dt = Omega must hold,
% nothing fancy.
%TT = real(Txpx.');
TT = real(TpBx.');

% Transform the input upsi to body coordinates and augment it with
% d2eta/dt2.
ixeta = Neta+[1:Neta].';
iueta = Ndj+[1:Neta].';  % Indices for d2eta/dt2 in u vector.
niueta = [[1:Ndj] Ndj+Neta+1:Nu].';
u = zeros(Nu,1);
u(niueta) = TpBu(niueta,niueta)*upsi;
u(iueta) = Txpx(ixeta,:)*dxpsi;

% Call the CLT function in body coordinates.
[L,R,y,A,B,C,D] =                                               ...
     buildClosedLoopTurbine (linFlag,x,u,s,a,epar,ppar,ypar,c,  ...
                             grav,P,shape0,mdamp0,              ...
                             Tas,Try,ch,Lel,foilwt,aoaz,aoast,  ...
                             xas,yas,Psi0,igen,ipit,iyaw);

W = xpsi(iW);

% (Think about using TT instead of real(TT) for term-by-term investigation.)

Lpsi = real(TT)*L*Txpx;  % real(TT) for comparing NL/Lin via complex step.
Rpsi = real(TT)*R;
ypsi = TBpy*y;

if (linFlag == 1)

   % Special consideration is needed when transforming d2eta/dt2.
   
   [dTpBx,dTBpx]   = derivMBC (Nx,blxdof(:,1),blxdof(:,2),blxdof(:,3),azi);
   [dTpBu,dTBpu]   = derivMBC (Nu,bludof(:,1),bludof(:,2),bludof(:,3),azi);
   [dTpBy,dTBpy]   = derivMBC (Ny,blydof(:,1),blydof(:,2),blydof(:,3),azi);

   Apsi = real(TT)*(A*Txpx - L*dTxpx + B(:,iueta)*dTxpx(ixeta,:));
   Apsi(:,iazi) = Apsi(:,iazi) ...
                + real(TT)*B(:,niueta)*dTpBu(niueta,niueta)*upsi;

   Bpsi = real(TT)*B*TpBu;
   Bpsi(:,iueta) = real(TT)*B(:,iueta)*Txpx(ixeta,ixeta);  % Overwrite OK.

   Cpsi = TBpy*(C*Txpx + D(:,iueta)*dTxpx(ixeta,:));
   Cpsi(:,iazi) = Cpsi(:,iazi) + dTBpy*y ...
                + TBpy*D(:,niueta)*dTpBu(niueta,niueta)*upsi;

   Dpsi = TBpy*D*TpBu;
   Dpsi(:,iueta) = TBpy*D(:,iueta)*Txpx(ixeta,ixeta);

else

   Apsi = sparse(size(A,1),size(A,2));
   Bpsi = sparse(size(B,1),size(B,2));
   Cpsi = sparse(size(C,1),size(C,2));
   Dpsi = sparse(size(D,1),size(D,2));

end
