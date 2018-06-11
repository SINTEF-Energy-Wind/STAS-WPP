function [dxadt,A,By,C,Dy,Psi] = aeroLin (psiFlag,xa,Vg,s,a,q,dqdt,P,bsh)
%
% Assemble the overall linear aerodynamic state space.  
%
%   States:           y vector:         u vector:
%   -------           ---------         ---------
%   ad        1       q         1:24       
%   a1,a2    2:3      dq/dt    25:48
%   Vih z,t  4:5      xng1,2   49:54
%   Vi z,t   6:7      vng1,2   55:60
%                     wg       61:63
%                     wa       64:66
%                     aq        67
%                     Cl,Cd,Cm 68:70
%                     Fl,Fd,M  71:73
%                     Fa       74:79
%                     Fp       80:85
%                     Fr       86:91
%                     Fzts     92:94
%                     Vg       95:97
%                     Ua       98:100
%                     Umag      101
%                     Ur      102:104
%                     Uzts    105:107
%                     Vzts    108:110
%                     Viq     111:112
%                     Viy     113:114
%                     Vixyz   115:117
%                     Wmag      118
%                     xeg     119:121
%                     xhg     122:124
%                     xnr1,2  125:130  
%                     xer     131:133
%                     r         124
%                     Lp        135
%                     z         136
%                     f         137
% (psiFlag = 1: repeat the above consecutively for each blade element
%               in an annulus.)
%                     Azi       (1)
%                     W         (1)
%                     Dp        (3)
%
% Version:        Changes:
% --------        -------------
% 20.01.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 20.01.2018      
%
% Inputs:
% -------
% psiFlag         : = 0 for independent calculation of dynamic Vi.
%                   = 1 for MBC calculation of dynamic Vi.
% x               : Vector of 7 aero states for each element.
% Vg              : 3*Nel, incoming wind, global coordinates.
% s,a             : Data structures.
% q, dqdt         : Structural DOFs.
% P               : Undeformed nodal positions.
% bsh             : If a.icp is negative, and blade modes are used,
%                   then this contains the mode shapes.
%
% Outputs:
% --------
% A,B,C,D         : State matrices.

Neb = a.Neb;
Nel = a.Nb*a.Neb;

Nx = 7*Nel;
Ny = 137*Nel + 5;
Nyae = 137;

% Allocate the full matrices.
Af  = spalloc (Nx,Nx,49*Nel);
Byf = spalloc (Nx,Ny,100*Nel);
Cf  = spalloc (Ny,Nx,50*Nel);
Dy  = spalloc (Ny,Ny,400*Nel);

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);

[iad,ia1,ia2,iVih,iVi,                    ...
 iq,idq,ixng1,ixng2,ivng1,ivng2,          ...
 iwg,iwa,iaq,iC,iFldm,iFa,iFp,Fr,iFzts,   ...
 iVg,iUa,iUmag,iUr,iUzts,iVzts,iViq,      ...
 iViy,iVixyz,iWmag,ixeg,ixhg,             ...
 ixnr1,ixnr2,ixer,irp,iLp,izPr,iPr] = getis (1);

% Compute transforms and other parameters needed for a BEM analysis.
[Tas,ch,Lel,foilwt,aoaz,aoast,xas,yas,iqel] = ...
                  BEMsetup (s,a,a.aoas,a.kfoils);
[Tar,Trg,TB0g,TsB,TBB0,dTar,dTsB,dTBB0,wg] = ...
                  BEMprepTransforms (s,a,q,dqdt,P,Tas);
[zr,Area,Dp,rp,Lp,xeg,xhg,xyg] = ...
                  BEMprepProjections (s,a,q,dqdt,P,Tn_y,Trg);
Dps = zeros(Nel,1);
Dps(1:Neb)         = Dp(1);
Dps(Neb+[1:Neb])   = Dp(2);
Dps(2*Neb+[1:Neb]) = Dp(3);

% Link wg with vng's.
for icomp = 1:3

   indr = iwg   + icomp + [0:Nyae:Nyae*(Nel-1)].';
   indc = ivng1 + icomp + [0:Nyae:Nyae*(Nel-1)].';
   Dy(indr,indc) = 0.5*eye(Nel);

   indr = iwg   + icomp + [0:Nyae:Nyae*(Nel-1)].';
   indc = ivng2 + icomp + [0:Nyae:Nyae*(Nel-1)].';
   Dy(indr,indc) = 0.5*eye(Nel);

end

% Extract the effective rotor speed from dqdt.  Note that this
% is stored in Dy elsewhere.
Ydof  = idofs(3);
Ddof  = idofs(4);
nodof = idofm(6) - 6;
[Omega,azi,jnk] = rotorSpeedAero (q,dqdt,P,Tn_y,Ydof,Ddof,nodof);

% Projected elements and diameter.  Indices indicate where the
% 45 y variables in projectElements lie in the sequence of 137 y
% variables associated with each element.
[y,Diajnk,ddype,ddDp] = projectElements (q,P,iqel,idofs,idofm,Tn_y);
% Rows: xeg, xnr1, xnr2, xer, r, Lp
indr = [ixeg+[1:3] ixnr1+[1:3] ixnr2+[1:3] ...
        ixer+[1:3] irp+1 iLp+1].';
% Cols: qe, xng1,2, xeg, xhg, xnr1,2, xer 
indc = [iq+[1:24] ixng1+[1:3] ixng2+[1:3] ixeg+[1:3] ...
        ixhg+[1:3] ixnr1+[1:3] ixnr2+[1:3] ixer+[1:3]].';
for iel = 1:Nel
   icN = Nyae*(iel-1);
   icr = size(indr,1)*(iel-1);
   icc = size(indc,1)*(iel-1);
   ir  = icN + indr;
   ic  = icN + indc;
   Dy(ir,ic) = Dy(ir,ic) + ddype(icr+[1:size(indr,1)],icc+[1:size(indc,1)]);
end
ir = [Ny-2:Ny].';
ic = [Nyae*(  Neb-1)+ixnr2+[1:3] ... % xnr2 of tip nodes.
      Nyae*(2*Neb-1)+ixnr2+[1:3] ...
      Nyae*(3*Neb-1)+ixnr2+[1:3]].';
Dy(ir,ic) = Dy(ir,ic) + ddDp;

% BEM.
Psi = aeroPsi (a,rp,bsh);
Psit = Psi.';
xaf = Psi*xa;
dxafdt = zeros(Nx,1);
for ian = 1:Neb

   val = Nyae;
   indN = [val*(ian-1)+[1:val]     ...
           val*(Neb+ian-1)+[1:val] ...
           val*(2*Neb+ian-1)+[1:val] [Ny-4:Ny]].';
   val = 72;
   ind72 = [val*(ian-1)+[1:val]     ...
            val*(Neb+ian-1)+[1:val] ...
            val*(2*Neb+ian-1)+[1:val]].';
   val = 45;
   ind45 = [val*(ian-1)+[1:val]     ...
            val*(Neb+ian-1)+[1:val] ...
            val*(2*Neb+ian-1)+[1:val]].';
   val = 36;
   ind36 = [val*(ian-1)+[1:val]     ...
            val*(Neb+ian-1)+[1:val] ...
            val*(2*Neb+ian-1)+[1:val]].';
   val = 9;
   ind9 = [val*(ian-1)+[1:val]     ...
           val*(Neb+ian-1)+[1:val] ...
           val*(2*Neb+ian-1)+[1:val]].';
   val = 7;
   ind7 = [val*(ian-1)+[1:val]     ...
           val*(Neb+ian-1)+[1:val] ...
           val*(2*Neb+ian-1)+[1:val]].';
   val = 3;
   ind3 = [val*(ian-1)+[1:val]     ...
           val*(Neb+ian-1)+[1:val] ...
           val*(2*Neb+ian-1)+[1:val]].';
   val = 2;
   ind2 = [val*(ian-1)+[1:val]     ...
           val*(Neb+ian-1)+[1:val] ...
           val*(2*Neb+ian-1)+[1:val]].';
   ind  = [ian Neb+ian 2*Neb+ian].'; 

   % Call BEMlin annulus-by-annulus.
   [dxafdt(ind7),aa,bby,cc,ddy] =                           ...
         BEMlin (psiFlag,xaf(ind7),                         ...
                 Tar(:,ind3),Tas(:,ind3),TsB(:,ind3),       ...
                 TBB0(:,ind3),TB0g(:,ind3),                 ...
                 dTar(:,ind72),dTsB(:,ind36),dTBB0(:,ind9), ...
                 Vg(ind3),wg(ind3),                         ...
                 zr(ind2),ch(ind),Lel(ind),                 ...
                 a.aoas,a.kfoils,foilwt(:,ind),             ...
                 aoaz(ind),aoast(ind2),                     ...
                 xas(ind),yas(ind),a.dens,                  ...
                 Dps(ind),rp(ind),Lp(ind),azi,Omega);

   % Store the annulus-by-annulus values in the overall matrices.
   Af(ind7,ind7)  = aa;
   Byf(ind7,indN) = bby;
   Cf(indN,ind7)  = cc;
   Dy(indN,indN) = Dy(indN,indN) + ddy;



end

dxadt = Psit*dxafdt;
A  = Psit*Af*Psi;
By = Psit*Byf;
C  = Cf*Psi;

