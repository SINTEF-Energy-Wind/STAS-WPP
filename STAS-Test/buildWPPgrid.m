function [dxdt,A,B] = buildWPPgrid (linFlag,x,u,gpar,xpar,nod,nsub)
%
%   States:              y vector:             u vector:
%                                              we          1
%                                              vpccd,q    2,3
%                                              itd,q     2*Nnod
% -------------------------- Collection grid -------------------------   
%   iLd,q  (2*Ncable)    we    (in) (u1)
%   vd,q   (2*Nnod)
%   ird,q  (2*Nnod)      ied,q (in) (2*Nnod)
%
% ----------------------------- AC export ----------------------------
%   i12d,q   1,2         we    (in) (u1)
%   ir2d,q   3,4         v1d,q (in) (u2,3)
%   v2d,q    5,6         v4d,q (in) (u4,5)
%   i23d,q   7,8
%   ir3d,q   9,10
%   v3d,q   11,12
%   i34d,q  13,14
%
% gpar            Ncab:  R     (Ohms)  Phase resistance in line j
%                 Ncab:  L     (H)     Phase inductance in line j
%                 Nnod:  C     (F)     Phase capacitance to ground at node k
% xpar            :  1:  a12   (-)     Np/Ns turns ratio.
%                  2-5:  I12   (A)     RMS phase currents for nonlinear inductance
%                  6-9:  L12   (H)     Lp + (a^2)Ls nonlinear phase inductances
%                   10:  R12   (Ohms)  Phase resistance
%                   11:  L2    (H)     Effective inductance,  Bus 2 -> ground
%                   12:  R2    (Ohms)  Effective resistance,  Bus 2 -> ground
%                   13:  C2    (F)     Effective capacitance, Bus 2 -> ground
%                   14:  L23   (H)     Phase inductance in cable
%                   15:  R23   (Ohms)  Phase resistance in cable
%                   16:  L3    (H)     Effective inductance,  Bus 3 -> ground
%                   17:  R3    (Ohms)  Effective resistance,  Bus 3 -> ground
%                   18:  C3    (F)     Effective capacitance, Bus 3 -> ground
%                   19:  a34   (-)     Np/Ns turns ratio.
%                20-23:  I34   (A)     RMS phase currents for nonlinear inductance
%                24-27:  L34   (H)     Lp + (a^2)Ls nonlinear phase inductances
%                   28:  R34   (Ohms)  Phase resistance
%

Ncab = size(nod,2);
Nnod = max(max(nod));

Nxg = 2*Ncab + 4*Nnod;  %number of states in collection grid
Nxx = 14;               %number of states in export system
Nxf  = Nxg + Nxx;

ixg = 0;
ixx = ixg + Nxg;

Nu  = 3 + 2*Nnod;

iut = 3;

% Identify locked DOFs at all the nodes, where there is no inductor path
% to ground. This is all the strings; only the export cable is compensated.
ldofs = [2*Ncab+2*Nnod+[1:2*Nnod] ixx+[3:4]].';
[ip,ret,cr] = partitionMatrix ([1:Nxf].',ldofs,[]);
Ndel = size(ldofs,1);
Nret = size(ret,1);
xf   = zeros(Nxf,1);
xf(ret) = x;

Lr = ones(Nnod,1);
Rr = zeros(Nnod,1);
gp = [gpar;Rr;Lr];
xg = xf(ixg+[1:Nxg]);

ug = [u(1);u(iut+[1:2*Nnod])];
isub = 1+2*(nsub-1)+[1:2];
ug(isub) = ug(isub) - xf(ixx+[1:2]);
[dxgdt,aag,bbg] = connectGrid (linFlag,xg,ug,gp,nod);

xx = xf(ixx+[1:Nxx]);
ux = [u(1);xf(ixg+2*Ncab+2*(nsub-1)+[1:2]);u(2:3)];
[dxxdt,aax,bbx] = ACexport (linFlag,xx,ux,xpar);

dxfdt = [dxgdt;dxxdt];
dxdt = dxfdt(ret);

if (linFlag == 1)

   Af = spalloc (Nxf,Nxf,0.1*Nxf*Nxf);
   Bf = spalloc (Nxf,Nu,0.1*Nxf*Nu);
   
   ir = ixg + [1:Nxg];
   ic = ir;
   Af(ir,ic) = aag;
   ic = 1;
   Bf(ir,ic) = Bf(ir,ic) + bbg(:,1);
   ic = ixx+[1:2];
   Af(ir,ic) = Af(ir,ic) - bbg(:,isub);
   ic = iut + [1:2*Nnod];
   Bf(ir,ic) = Bf(ir,ic) + bbg(:,1+[1:2*Nnod]);

   ir = ixx + [1:Nxx];
   ic = ir;
   Af(ir,ic) = aax;
   ic = 1;
   Bf(ir,ic) = Bf(ir,ic) + bbx(:,1);
   ic = ixg+2*Ncab+2*(nsub-1)+[1:2];
   Af(ir,ic) = Af(ir,ic) + bbx(:,2:3);
   ic = [2:3];
   Bf(ir,ic) = Bf(ir,ic) + bbx(:,4:5);

else

   Af = sparse(Nxf,Nxf);
   Bf = sparse(Nxf,Nu);

end

A = Af(ret,ret);
B = Bf(ret,:);
