function [L,R,yout,A,B,C,D] =                                              ...
              buildClosedLoopTurbine (linFlag,x,u,s,a,epar,ppar,ypar,c,    ...                                     
                                      grav,P,shape0,mdamp0,                ...
                                      Tas,Try,ch,Lel,foilwt,aoaz,aoast,    ...
                                      xas,yas,Psi0,igen,ipit,iyaw)
%
% Link the turbine control functions with the wind turbine.
%
%   States:              y vector:                u vector:
%                                                 F       Ndj   (Env)
%                                              d2eta/dt2  Neta  (u)
%                                                 Vg    3*Nae   (Env)
%                                                 we        1   (Grid)
%                                                 th_e      2   (Grid)
%                                                 vsd,q   3,4   (Grid)
%                                                 Pc        5   (Pl.cont.)
%                                                 Qc        6   (Pl.cont.)
% ----------------------- Structure --------------------------
%   eta        N         q         Ndj            F       Ndj   (u)
%   deta/dt    N         dq/dt     Ndj         d2eta/dt2  Neta  (u)
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
%   ba         1         b,dbdt    (-)  (Turb)    bhat      1   (Cont)
%   dbadt      2         Ta         1
%
% ----------------------- Electrical --------------------------
%   igd,q     1,2        wg        (-)  (Turb)    we        1   (u)
%   Vdc        3         Tg         1             th_e      2   (u)
%   ipd,q     4,5        isd,q     2,3            vsd,q    3,4  (u)
%   imgd,q    6,7                                 ihgd,q   5,6  (Cont)
%   Psig      8,9                                 Vhdc      7   (cpar)
%   wemg      10                                  Qh        8   (u)
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
% 02.02.2019      Original code.
% 23.01.2020      Minor updates to accommodate the revisions to the
%                 structural equations-of-motion and the u vector.
%
% Version:        Verification:
% --------        -------------
% 02.02.2019      
% 23.01.2020      Derivatives of A matrix verified via complex step.
%
% Inputs:
% -------
% 
%
% Outputs:
% --------
% 

% Determine some indices.
[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
Nnod = Ndof/6;
Ndj  = Ndof + 6;
Nae  = a.Nb*a.Neb;
Neta = size(shape0,2);
Nxs  = 2*Neta;
Nxa  = size(Psi0,2);

Nxp = 6;
Nxy = 2;
Nxe = 25;
Nxc = 31;

Nxt = Nxs + Nxa + Nxp + Nxy + Nxe;
Nx1  = 2*Neta + Nxa;

Nys  = 4*Ndj + 9*Nnod;
Nya  = 137*Nae + 5;
Ny1  = Nys + Nya;
Ny1n = 4*Ndj;     % The nonlinear outputs.

% Note, the indexing here is for the u vector input.  Note the
% difference in Pc, Qc (u vector) versus Qh, Pc (matrices).
Nuin = size(u,1);
Nus  = Ndj + Neta;
Nua  = 3*Nae;
Nu1  = Nus + Nua;
Nup  = 0;
Nue  = 4;
Nuc  = 2;

ius  = 0;
iua  = ius + Nus;
iup  = iua + Nua;
iue  = iup + Nup;
iuc  = iue + Nue;

F    = u(ius+[1:Ndj]);
etadd= u(ius+Ndj+[1:Neta]);
Vg   = u(iua+[1:3*Nae]);
we   = u(iue+1);
th_e = u(iue+2);
vs   = [u(iue+3);u(iue+4)];
Pc   = u(iuc+1);
Qc   = u(iuc+2);

% Call the nonlinear control first.  The controller command outputs are
% functions of the measured rotor azimuth (for IBP), and the internal 
% controller states.  So provided that the azimuth can be determined,
% the control commands will be correct, regardless of whether the other
% inputs are in error.  For the azimuth, use the driveshaft body 
% rotary joint position (rear bearing).  This may differ subtly from
% the location where the rotor speed is measured, but it is expedient
% to use this definition of azimuth here.
azi   = x(Neta-4);
jnkW  = x(2*Neta-4);
bet   = x(Neta-[2:-1:0]);
jnkPe = Pc;
jnkvT = [0;0];
jnkMbl = [0;0;0];
jnkwang = 0;
xc = x(Nx1+Nxp+Nxy+Nxe+[1:Nxc]);
uc = [Pc;azi;jnkW;bet;jnkPe;jnkvT;jnkMbl;jnkwang];

[jnkdx,ycout,jnka,jnkb,jnkc,jnkd,blyc,bluc] =                        ...
        turbineControl (0,xc,uc,c.cpar,c.cpct,                       ...
                        c.KeTab,c.WVTab,c.WPTab,c.bminTab,c.KTables, ...
                        c.KFTab,c.KSTab,c.KSqTab,c.KpiTab,c.KiiTab,c.RSCFlag);

% Now we have the correct control commands.  Call the linear/nonlinear
% turbine model.
bhat = ycout(1:3);
ihg  = ycout(4:5);
yang = 0;    % Yaw control not yet implemented.
xt = x(1:Nx1+Nxp+Nxy+Nxe);
ut = [F;etadd;Vg;bhat;yang;we;th_e;vs;ihg;c.Vhdc;Qc];
[llt,rrt,ytout,aat,bbt,cct,ddt] =                                 ...
      buildOpenLoopTurbine (linFlag,xt,ut,s,a,epar,ppar,ypar,     ...
                            grav,P,shape0,mdamp0,                 ...
                            Tas,Try,ch,Lel,foilwt,aoaz,aoast,     ...
                            xas,yas,Psi0,igen,ipit,iyaw);

% Call the linear/nonlinear controller again to get the correct dx/dt
% and if needed the linearized state matrices.
q  = ytout(1:Ndj);
dq = ytout(Ndj+[1:Ndj]);
[Wg,ddyg] = genSpeed (q(igen),dq(igen),P(igen));
Pe = (ytout(Ny1n+4+[2:3]).')*u(Nu1+[3:4]);
TB0g = TFromTheta (P(idofs(3)+[4:6]));
[TBB0,dTBB0] = dTdth (q(idofs(3)+[4:6]));
TgB = (TB0g*TBB0).';
vT = TgB(1:2,:)*dq(idofs(3)+[1:3]) + dq(idofs(3)+6+[1:2]);
wang = 0;    % Yaw control not yet implemented.

zetpy = zeros(3,1);  % Use nodal rotations as a proxy for blade root moments.
for ib = 1:3
   qn1 = q(idofs(5+ib)+6+[1:6]);
   Pn1 = P(idofs(5+ib)+6+[1:6]);
   qn2 = q(idofs(5+ib)+12+[1:6]);
   Pn2 = P(idofs(5+ib)+12+[1:6]);
   [xe,TsB] = elementCSFromNodes (qn1,qn2,Pn1,Pn2);
   mu = getMu (qn1,qn2,Pn1,Pn2,TsB);
   zets = mu(4:6);
   zetpy(ib) = TsB(2,:)*zets;
end

uc = [Pc;azi;Wg;bet;Pe;vT;zetpy;0];

[dxcdt,ycout,aac,bbc,ccc,ddc,blyc,bluc] =                            ...
        turbineControl (linFlag,xc,uc,c.cpar,c.cpct,                 ...
                        c.KeTab,c.WVTab,c.WPTab,c.bminTab,c.KTables, ...
                        c.KFTab,c.KSTab,c.KSqTab,c.KpiTab,c.KiiTab,c.RSCFlag);

L = [llt sparse(Nxt,Nxc);sparse(Nxc,Nxt) speye(Nxc)];
R = [rrt;dxcdt];
yout = [ytout;ycout];

Nx = size(R,1);
Nu = size(u,1);
Ny = size(cct,1) + size(ccc,1);

if (linFlag == 1)

   % Define a z vector for linking, consisting of the local u variables,
   % and then a subset of the local y variables: q and dq/dt, systems,
   % and controls.
   ixs = 0;
   ixa = ixs + Nxs;
   ixp = ixa + Nxa;
   ixy = ixp + Nxp;
   ixe = ixy + Nxy;
   ixc = ixe + Nxe;

   Nus = Ndj + Neta;
   Nua = 3*Nae;
   Nup = 3;
   Nuy = 1;
   Nue = 8;
   Nuc = 13;
   Nut = Nus + Nua + Nup + Nuy + Nue;
   ius = 0;
   iua = ius + Nus;
   iup = iua + Nua;
   iuy = iup + Nup;
   iue = iuy + Nuy;
   iuc = iue + Nue;

   Nys = 4*Ndj;
   Nya = 0;
   Nyp = 3;
   Nyy = 1;
   Nye = 3;
   Nyc = 6;
   Nyt = Nys + Nya + Nyp + Nyy + Nye;
   iys = iuc + Nuc;
   iyp = iys + Nys;
   iyy = iyp + Nyp;
   iye = iyy + Nyy;
   iyc = iye + Nye;

   % (This should be done in a cleaner way... I've dropped some y's in order
   %  that the linear matrix row entries match the NL vector entries...)
   jz1 = [[1:3*Ndj] 3*Ndj+9*Nnod+[1:Ndj] Ny1+[3 6 9] Ny1+9+3 Ny1+12+[2 3 4]];

   Nz = iyc + Nyc;

   Az = spalloc (Nx,Nx,0.3*Nx*Nx);
   Bu = spalloc (Nx,Nuin,0.2*Nx*Nuin);
   Bz = spalloc (Nx,Nz,0.1*Nx*Nz);
   Cz = spalloc (Nz,Nx,0.1*Nz*Nx);
   Du = spalloc (Nz,Nuin,0.05*Ny*Nuin);
   Dz = spalloc (Nz,Nz,0.02*Nz*Nz);

   Az(ixs+[1:Nxt],ixs+[1:Nxt]) = aat;
   Bz(ixs+[1:Nxt],ius+[1:Nut]) = bbt;
   Cz(iys+[1:Nyt],ixs+[1:Nxt]) = cct(jz1,:);
   Dz(iys+[1:Nyt],ius+[1:Nut]) = ddt(jz1,:);

   Az(ixc+[1:Nxc],ixc+[1:Nxc]) = aac;
   Bz(ixc+[1:Nxc],iuc+[1:Nuc]) = bbc;
   Cz(iyc+[1:Nyc],ixc+[1:Nxc]) = ccc(1:Nyc,:);
   Dz(iyc+[1:Nyc],iuc+[1:Nuc]) = ddc(1:Nyc,:);

   % Link z vector u's with global u's.
   ir = [ius+[1:Nus] iua+[1:Nua] iue+[1 2 3 4 8] iuc+1];
   ic = [[1:Ndj+Neta] Ndj+Neta+[1:3*Nae] Nu1+[1 2 3 4 6 5]];
   ss = ones(1,Ndj+Neta+3*Nae+6);
   Du = Du + sparse (ir,ic,ss,Nz,Nuin);

   % Link z vector system u's with z vector control y's.
   ir = [iup+[1 2 3] iuy+1 iue+[5 6]];
   ic = [iyc+[1 2 3 6] iyc+[4 5]];
   ss = ones(1,6);
   Dz = Dz + sparse (ir,ic,ss,Nz,Nz);

   % Link the (simple) z vector control u's with the turbine x's.
   ir = iuc+[2 4 5 6];
   ic = [ixs+Neta-[4 2 1 0]];
   ss = ones(1,4);
   Cz = Cz + sparse(ir,ic,ss,Nz,Nx);

   % Linearization of rotor speed sensor input to control.
   ir = iuc + 3;
   ic = iys+[igen Ndj+igen];
   Dz(ir,ic) = Dz(ir,ic) + ddyg;

   % Linearization of electric power.  Pe = isd*vsd + isq*vsq.
   ir = iuc + [7 7];
   ic = Nu1 + [3 4];
   ss = yout(iye-iys+[2 3]);
   Du = Du + sparse (ir,ic,ss,Nz,Nuin);
   ic = iye + [2 3];
   ss = u(Nu1+[3 4]);
   Dz = Dz + sparse (ir,ic,ss,Nz,Nz);

   % Linearization of vT.
   % vT = TgB(1:2,:)*dq(idofs(3)+[1:3]) + dq(idofs(3)+6+[1:2]);
   ir = iuc + [8 9];
   ic = iys + Ndj + idofs(3) + 6 + [1 2];
   Dz(ir,ic) = Dz(ir,ic) + speye(2);
   ic = iys + Ndj + idofs(3) + [1:3];
   Dz(ir,ic) = Dz(ir,ic) + TgB(1:2,:);
   for jj = 1:3
      jc3 = 3*(jj-1);
      ic = iys + idofs(3) + 3 + jj;
      TT = (TB0g*dTBB0(:,jc3+[1:3])).';
      Dz(ir,ic) = Dz(ir,ic) + TT(1:2,:)*dq(idofs(3)+[1:3]);
   end

   % Linearization of blade root rotations.
   for ib = 1:3
      qn1 = q(idofs(5+ib)+6+[1:6]);
      Pn1 = P(idofs(5+ib)+6+[1:6]);
      qn2 = q(idofs(5+ib)+12+[1:6]);
      Pn2 = P(idofs(5+ib)+12+[1:6]);
      [xe,TsB] = elementCSFromNodes (qn1,qn2,Pn1,Pn2);
      mu = getMu (qn1,qn2,Pn1,Pn2,TsB);
      dTsB = derivElementCS (qn1,qn2,Pn1,Pn2,TsB);
      dmu = dmudq (qn1,qn2,Pn1,Pn2,TsB,dTsB);
      ir = iuc + 9 + ib;
      ic = idofs(5+ib) + 6 + [1:12];
      Dz(ir,ic) = Dz(ir,ic) + TsB(2,:)*dmu(4:6,:);
      for jj = 1:12
         jc3 = 3*(jj-1);
         Dz(ir,ic(jj)) = Dz(ir,ic(jj)) + dTsB(2,jc3+[1:3])*mu(4:6);
      end
   end

   % Reduce the matrices.
   Dnorm = ones(Nz,1);
   [A,B,C,D] = modularToUnifiedStateSpace (Az,Bu,Bz,Cz,Du,Dz,Dnorm);
   C = C(iys+1:Nz,:);
   D = D(iys+1:Nz,:);

else

   Nyo = size(yout,1);
   A = sparse (Nx,Nx);
   B = sparse (Nx,Nuin);
   C = sparse (Nyo,Nx);
   D = sparse (Nyo,Nuin);

end

