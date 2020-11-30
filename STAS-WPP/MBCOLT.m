function [Lop,Rop,yop,Aop,Bop,Cop,Dop] =                    ...
                  MBCOLT (linFlag,xop,dxop,uop,s,a,         ...
                          epar,bpar,ypar,m,c,               ...
                          grav,P,shape0,mdamp0,             ...
                          Tas,Try,ch,Lel,foilwt,aoaz,aoast, ...
                          xas,yas,Psi0,igen,ipit,iyaw)
%
% Return the multi-blade coordinate transformed version of the dynamic
% open-loop turbine equations.  This adds sensors and a generator
% power control function to the basic buildOpenLoopTurbine.
%
% Caution, the ighq input is unused, it is derived within based on the
% generator power control.
%
%   States:              y vector:             u vector:
% ----------------------- Structure --------------------------
%   eta        N         q         Ndj            F       Ndj   (Env)
%   deta/dt    N         dq/dt     Ndj         d2eta/dt2  Neta  (ext. solution)
%                        d2q/dt2   Ndj
%                        xng     3*Nnod
%                        vng     6*Nnod
%                        F         Ndj
%
% ------------208------- Aerodynamic -------------------------
%   ad         1         (q)       1:24           Vg        3   (Env)
%   a1,a2     2:3        (dq/dt)  25:48
%   Vih z,t   4:5        xng1,2   49:54
%   Vi z,t    6:7        vng1,2   55:60
%                        wg       61:63
%                        wa       64:66
%                        aq        67
%                        Cl,Cd,Cm 68:70
%                        Fl,Fd,M  71:73
%                        Fa       74:79
%                        Fp       80:85
%                        Fr       86:91
%                        Fzts     92:94
%                        Vg       95:97
%                        Ua       98:100
%                        Umag      101
%                        Ur      102:104
%                        Uzts    105:107
%                        Vzts    108:110
%                        Viq     111:112
%                        Viy     113:114
%                        Vixyz   115:117
%                        Wmag      118
%                        xeg     119:121
%                        xhg     122:124
%                        xnr1,2  125:130  
%                        xer     131:133
%                        r         134
%                        Lp        135
%                        z         136
%                        f         137  
% (Repeat the above x,y,u consecutively for each blade element.)
%                        Azi       (1)  (effective rotor azimuth)
%                        Waero     (1)  (aero rotor speed)
%                        Dp        (3)
%
% ------------334--------- Actuators --------------------------
% (Repeat for each blade and yaw system.)
%   ba         1         b,dbdt    1,2  (Turb)     bhat      1   (Cont)
%   dbadt      2         Ta         3
%
% ------------342-------- Electrical --------------------------
%   igd,q     1,2        wg         1   (Turb)     we        1   (Grid)
%   Vdc        3         Tg         2              th_e      2   (Grid)
%   ipd,q     4,5        isd,q     3,4             vsd,q    3,4  (Grid)
%   imgd,q    6,7                                  ihgd,q   5,6  (Cont)
%   Psig      8,9                                  Vhdc      7   (Cont)
%   wemg      10                                   Qh        8   (Cont)
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
% ------------367---------- Sensors ---------------------------
%   Wm         1         W          1
%   bm        2:4        bet       2:4
%   ym         5         yaw        5
%   Pem        6         Pe         6
%   vm        7,8        vnac      7,8
%   Vam        9         Va         9
%   thm       10         tha       10
%
% ------------377--------- P control --------------------------
%  (Pem)                 ighq       1              Phat      1
%   PsiP       1         Pe         2             (Pe)
%   nd        2,3
%
% Version:        Changes:
% --------        -------------
% 12.02.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 12.02.2020      
%

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
Ndj  = Ndof + 6;
Nnod = Ndof/6;
[imdofs,Neta] = getmdofRefs (s);
Neb  = s.blade(1).Nel;
Nae  = a.Nb*Neb;

Nxs  = 2*Neta;
Nus  = Ndj + Neta;
Nys  = 4*Ndj;
Nxa  = size(Psi0,2);
Nua  = 3*Nae;
Nya  = 0;
Nxb  = 6;
Nub  = 3;
Nyb  = 3;
Nxy  = 2;
Nuy  = 1;
Nyy  = 1;
Nxe  = 25;
Nue  = 8;
Nye  = 3;
Nxm  = 10;
Num  = 0;
Nym  = 10;
Nxc  = 3;
Nuc  = 1;
Nyc  = 2;
Nysl = 4*Ndj + 9*Nnod; % Additional y indices in bOLT matrices.
Nyal = 137*Nae + 5;
Nybl = 9;
Nyyl = 3;
Nyel = 4;

jxs  = 0;
jus  = 0;
jys  = 0;
jxa  = jxs + Nxs;
jua  = jus + Nus;
jya  = jys + Nys;
jxb  = jxa + Nxa;
jub  = jua + Nua;
jyb  = jya + Nya;
jxy  = jxb + Nxb;
juy  = jub + Nub;
jyy  = jyb + Nyb;
jxe  = jxy + Nxy;
jue  = juy + Nuy;
jye  = jyy + Nyy;
jxm  = jxe + Nxe;
jum  = jue + Nue;
jym  = jye + Nye;
jxc  = jxm + Nxm;
juc  = jum + Num;
jyc  = jym + Nym;

jyal = jys  + Nysl;
jybl = jyal + Nyal;
jyyl = jybl + Nybl;
jyel = jyyl + Nyyl;
jyml = jyel + Nyel;
jycl = jyml + Nym;

iynl = [jys+[1:3*Ndj] jys+3*Ndj+9*Nnod+[1:Ndj] jybl+[3 6 9 12] ...
        jyel+[2 3 4] jyml+[1:Nym] jycl+[1:Nyc]].';

Nx   = size(xop,1);
Nu   = size(uop,1);
Ny   = jyc + 2;
Nyl  = jycl + 2;

iW   = 2*Neta - 4;
iazi = Neta - 4;
azi  = xop(iazi);

Nxa3  = Nxa/3;
Nyal3 = 137*Neb;
Neb3  = 3*Neb;

[b1eta,b2eta,b3eta] = MBCindices_Nmd (Neta,imdofs);
[b1q,b2q,b3q]       = MBCindices_Ndj (Ndj,idofs);
bn3 = [3*(inods(6)-1)+[1:3*s.blade(1).Nnod].' ...
       3*(inods(7)-1)+[1:3*s.blade(1).Nnod].' ...
       3*(inods(8)-1)+[1:3*s.blade(1).Nnod].'];
bn6 = [6*(inods(6)-1)+[1:6*s.blade(1).Nnod].' ...
       6*(inods(7)-1)+[1:6*s.blade(1).Nnod].' ...
       6*(inods(8)-1)+[1:6*s.blade(1).Nnod].'];
baero = [[1:Nyal3].' Nyal3+[1:Nyal3].' 2*Nyal3+[1:Nyal3].'];

blx = [jxs+[b1eta b2eta b3eta];                            ...
       jxs+Neta+[b1eta b2eta b3eta];                       ...
       jxa+[[1:Nxa3].' Nxa3+[1:Nxa3].' 2*Nxa3+[1:Nxa3].']; ...
       jxb+[[1;2] [3;4] [5;6]];                            ...
       jxm+[2 3 4]];

bly = [jys+[b1q b2q b3q];                     ...
       jys+Ndj+[b1q b2q b3q];                 ...
       jys+2*Ndj+[b1q b2q b3q];               ...
       jys+3*Ndj+bn3;                         ...
       jys+3*Ndj+3*Nnod+bn6;                  ...
       jys+3*Ndj+9*Nnod+[b1q b2q b3q];        ...
       jyal+baero;                            ...
       jybl-[2 1 0];                          ...  % Dp's.
       jybl+[[1 2 3].' [4 5 6].' [7 8 9].'];  ...
       jyml+[2 3 4]];

blu = [jus+[b1q b2q b3q];                                ...
       jus+Ndj+[b1eta b2eta b3eta];                      ...
       jua+[[1:Neb3].' Neb3+[1:Neb3].' 2*Neb3+[1:Neb3].']; ...
       jub+[1 2 3]];

[TpBx,TBpx]     = MBC (Nx,blx(:,1),blx(:,2),blx(:,3),azi);
[TpBu,TBpu]     = MBC (Nu,blu(:,1),blu(:,2),blu(:,3),azi);
[TpBy,TBpy]     = MBC (Nyl,bly(:,1),bly(:,2),bly(:,3),azi);
[dTpBx,dTBpx]   = derivMBC (Nx,blx(:,1),blx(:,2),blx(:,3),azi);

Txpx = TpBx;
Txpx(:,iazi) = Txpx(:,iazi) + dTpBx*xop;

TT = TpBx.';

xo = TpBx*xop;

iueta = Ndj+[1:Neta].';  % Indices for d2eta/dt2 in u vector.
niueta = [[1:Ndj] Ndj+Neta+1:Nu].';
ixeta = Neta+[1:Neta].';
uo = zeros(Nu,1);
uo(niueta) = TpBu(niueta,niueta)*uop(niueta);
uo(iueta) = Txpx(ixeta,:)*dxop;

% Do an initial solution of genPcontrol without knowing the iis*vvs
% power.  This will allow us to obtain the current command to input
% to buildOpenLoopTurbine.  It is OK not to know the input Pe, since
% the output is not a direct function of this.
nnxc = 3;  % 4 states in genPcontrol... Pem is assigned to the sensor block.
nnyc = 2;  % Going to shuffle Pe from u to y...
nnuc = 1;
indx = [jxm+6, jxc+[1:nnxc]].';
gpar = [c.Kp;c.Ki;m.aP;c.anp;c.z1np;c.z2np];
xc   = xop(indx);
uc   = [uop(juc+1); 0];
[jnkdxcdt,ighq,jnkaac,jnkbbc,jnkccc,jnkddc] = genPcontrol (xc,uc,gpar);

% buildOpenLoopTurbine does not include the sensors or P control.
% Account for P control on the ighq input.
xin = xo(1:jxm);
uin = uo(1:jum);
uin(jue+6) = ighq;

[Lo1,Ro1,yo1,Ao1,Bo1,Co1,Do1] =                                 ...
      buildOpenLoopTurbine (linFlag,xin,uin,s,a,epar,bpar,ypar, ...
                            grav,P,shape0,mdamp0,               ...
                            Tas,Try,ch,Lel,foilwt,aoaz,aoast,   ...
                            xas,yas,Psi0,igen,ipit,iyaw);

% Prepare the generator power control, so that we are dealing with
% an electrical power command, rather than a converter current
% command.  Use the LP filter parameter from the sensor block,
% rather from the old control block.
iis  = yo1(jye+[2 3]);
vvs  = uop(jue+[3 4]);
Pes  = iis(1)*vvs(1) + iis(2)*vvs(2);
xc   = xop(indx);
uc   = [uop(juc+1); Pes];
[dxcdt,ycout,aac,bbc,ccc,ddc] = genPcontrol (xc,uc,gpar);

% Prepare the sensors (measurements, denoted "m").
mvec = [m.aW; m.ab; m.ab; m.ab; m.ay; m.aP; m.ad; m.ad; m.aV; m.at];
nnxm = size(mvec,1);
nnym = size(mvec,1);
nnum = 0;

% Note, the y vector is based on the "reduced" indexing.
Lo = sparse(Nx,Nx);
Ro = zeros(Nx,1);
yo = zeros(Ny,1);
Lo(1:jxm,1:jxm) = Lo1;
Ro(1:jxm) = Ro1;
yo(1:jym) = yo1;

% Link in the sensors.  Let the nacelle accleration (velocity)
% be measured at the tower master node, in global coordinates.
% The y-vector measured power should come via the generator P
% control, not the sensor Pe, which should be zero.
qB  = yo1(idofs(2)+[1:6]);
qBd = yo1(Ndj+idofs(2)+[1:6]);
PB  = P(idofs(2)+[1:6]);
qn  = yo1(idofm(3)+[1:6]);
qnd = yo1(Ndj+idofm(3)+[1:6]);
Pn  = P(idofm(3)+[1:6]);
[vnac,ddyv] = globalVelocity (qn,qB,Pn,PB,qnd,qBd);
Vwnd = sqrt(uop(jua+1)^2+uop(jua+2)^2);
thw = atan2c(uop(jua+2),uop(jua+1));
Lo(jxm+[1:nnxm],jxm+[1:nnxm]) = speye(nnxm);
rvec = [xo(iW); xo(jxs+Neta-[2 1 0]); xo(jxs+Neta-5);  ...
        Pes; vnac(1:2); Vwnd; thw];
Ro(jxm+[1:nnxm]) = -mvec.*xo(jxm+[1:nnxm]) + mvec.*rvec;
yo(jym+[1:nnym]) = rvec;

% Link in the generator P control, overwriting the existing sensor
% LP filter.
Lo(jxc+[1:nnxc],jxc+[1:nnxc]) = speye(nnxc);
Ro(indx) = dxcdt;
yo(jyc+1) = ycout;
yo(jyc+2) = Pes;
yo(jym+6) = 0;

Lop = real(TT)*Lo*Txpx;
Rop = real(TT)*Ro;
yop = TBpy(iynl,iynl)*yo;

if (linFlag == 1)

   % Note, these intermediate matrices are based on the "full" y
   % vector indexing.
   AA  = sparse(Nx,Nx);
   BBu = sparse(Nx,Nu);
   BBy = sparse(Nx,Nyl);
   CC  = sparse(Nyl,Nx);
   DDu = sparse(Nyl,Nu);
   DDy = sparse(Nyl,Nyl);
   Lop(1:jxm ,1:jxm) = Lo1;
    AA(1:jxm ,1:jxm) = Ao1;
   BBu(1:jxm ,1:jum) = Bo1;
    CC(1:jyml,1:jxm) = Co1;
   DDu(1:jyml,1:jum) = Do1;

   % Sensors.  Do not include measured power, this will come from
   % the generator P control.
   inds = [[1:5] [7:10]].';
   aam = -diag(mvec);
   bbym = diag(mvec);
   AA(jxm+inds,jxm+inds) = aam(inds,inds);  % Sensor dynamics.
   BBy(jxm+inds,jyml+inds) = bbym(inds,inds); % Effect of sensor y on sensor x.
   CC(jyml+1,iW) = 1;                        % Link driveshaft W to sensor W.
   CC(jyml+[2 3 4],jxs+Neta-[2 1 0]) = speye(3); % Link pitch beta to sensor beta.
   CC(jyml+5,jxs+Neta-5) = 1;                % Link yaw angle to sensor yaw.
%   DDy(jyml+6,jyel+[3 4]) = [vvs(1), vvs(2)];  % Link is to sensor Pe.
%   DDu(jyml+6,jue+[3 4]) = [iis(1), iis(2)]; % Link vs to sensor Pe.
   DDy(jyml+[7 8],idofs(2)+[1:6]) = ddyv(1:2,1:6); % Link qB to sensor vnac.
   DDy(jyml+[7 8],idofm(3)+[1:6]) = ddyv(1:2,7:12); % Link qn to sensor vnac.
   DDy(jyml+[7 8],Ndj+idofs(2)+[1:6]) = ddyv(1:2,13:18); % Link qBd to sensor vnac.
   DDy(jyml+[7 8],Ndj+idofm(3)+[1:6]) = ddyv(1:2,19:24); % Link qnd to sensor vnac.
   DDu(jyml+9,jua+[1 3*Neb+1 6*Neb+1]) = uop(jua+1)/(3*Vwnd); % Link Vg0x at inner blade nodes to sensor Va.
   DDu(jyml+9,jua+[2 3*Neb+2 6*Neb+2]) = uop(jua+2)/(3*Vwnd); % Link Vg0y at inner blade nodes to sensor Va.
   DDu(jyml+10,jua+[1 3*Neb+1 6*Neb+1]) =-uop(jua+2)/(3*Vwnd^2); % Link Vg0x at inner blade nodes to sensor thw.
   DDu(jyml+10,jua+[2 3*Neb+2 6*Neb+2]) = uop(jua+1)/(3*Vwnd^2); % Link Vg0y at inner blade nodes to sensor thw.

   % Link in the generator P control.
   indx = [jxm+6, jxc+[1:3]].';
   AA(indx,indx) = aac;                      % Gen P control dynamics.
   BBu(indx,juc+1) = bbc(:,1);               % Link Phat to P control states.
   CC(jycl+1,indx) = ccc;                    % Link states to P control ighq.
   DDu(jycl+1,juc+1) = ddc(:,1);             % Link Phat to P control ighq.

   BBy(indx,jycl+2) = bbc(:,2);              % Link Pe to control states.
   DDy(jycl+1,jycl+2) = ddc(:,2);            % Link Pe to P control ighq.
   BBy(:,jycl+1) = BBu(:,jue+6);             % P control ighq behaves as elec ighq.
   DDy(:,jycl+1) = DDu(:,jue+6);
   DDy(jycl+2,jyel+[3 4]) = [vvs(1), vvs(2)]; % Link is to P control Pe.
   DDu(jycl+2,jue+[3 4]) = [iis(1), iis(2)]; % Link vs to P control Pe.

   % Unhook the ighq input, to be consistent with the nonlinear part,
   % where this has been overridden by the P control.
   BBu(:,jue+6) = sparse(Nx,1);
   DDu(:,jue+6) = sparse(Nyl,1);

   [Ao,Bo,Co,Do] = modularToUnifiedStateSpace (AA,BBu,BBy,CC,DDu,DDy,ones(Nyl,1));
   Nxop = size(Ao,1);
   Nuop = size(Bo,2);
   Nyop = size(Co,1);

   % Linearized MBC transforms.
   [dTpBu,dTBpu]   = derivMBC (Nu,blu(:,1),blu(:,2),blu(:,3),azi);
   [dTpBy,dTBpy]   = derivMBC (Nyl,bly(:,1),bly(:,2),bly(:,3),azi);
   [d2TpBx,d2TBpx] = secondDerivMBC (Nx,blx(:,1),blx(:,2),blx(:,3),azi);

   W = xop(iW);
   dTxpx = W*dTpBx;
   dTxpx(:,iazi) = dTxpx(:,iazi) + dTpBx*dxop + W*d2TpBx*xop;

   Aop = real(TT)*(Ao*Txpx - Lo*dTxpx + Bo(:,iueta)*dTxpx(ixeta,:));
   Aop(:,iazi) = Aop(:,iazi) ...
               + real(TT)*Bo(:,niueta)*dTpBu(niueta,niueta)*uop(niueta);
   Bop = real(TT)*Bo*TpBu;
   Bop(:,iueta) = real(TT)*Bo(:,iueta)*Txpx(ixeta,ixeta);  % Overwrite OK.

   Cop = TBpy(iynl,:)*(Co*Txpx + Do(:,iueta)*dTxpx(ixeta,:));
   Cop(:,iazi) = Cop(:,iazi) + dTBpy(iynl,iynl)*yo ...
                + TBpy(iynl,:)*Do(:,niueta)*dTpBu(niueta,niueta)*uop(niueta);

   Dop = TBpy(iynl,:)*Do*TpBu;
   Dop(:,iueta) = TBpy(iynl,:)*Do(:,iueta)*Txpx(ixeta,ixeta);

else

   Aop = sparse (Nx,Nx);
   Bop = sparse (Nx,Nu);
   Cop = sparse (Ny,Nx);
   Dop = sparse (Ny,Nu);

end

