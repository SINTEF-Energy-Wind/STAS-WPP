function [dxdt,A,By,C,Dy] = BEMlin (psiFlag,x,                     ...
                                    Tar,Tas,TsB,TBB0,TB0g,         ...
                                    dTar,dTsB,dTBB0,               ...
                                    Vg,wg,zr,ch,Lel,               ...
                                    aoas,kfoils,foilwt,aoaz,aoast, ...
                                    xas,yas,rho,Dp,rp,Lp,azi,Omega)
%
% Assemble the linear aerodynamic state space.  This could be called
% with all elements at once, but it's probably better, in light of
% the large y vector, to call it element-by-element, or if psiFlag
% = 1, annulus-by-annulus.
%
% dx/dt is included as an output for debugging by comparing with
% BEMNL.m.  For nonlinear time domain simulation, use BEMNL.m, as
% it is much faster.
%
% The structural DOFs q consist of qy, qB, qn1, and qn2.  
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
%                     r         134
%                     Lp        135
%                     z         136
%                     f         137
% (Repeat the above consecutively for each blade element.)
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
% 20.01.2018      Verified that the dx/dt output matches that of BEMNL.
%                 Verified that the unified A matrix, A + By inv(I-Dy) C,
%                 matches the complex step values for x input, dx/dt output.
%                 [Need to verify annulus indexing when psiFlag = 1. Need
%                  to verify Fp.]
%
% Inputs:
% -------
% psiFlag         : = 0 for independent calculation of dynamic Vi.
%                   = 1 for MBC calculation of dynamic Vi.
% x               : Vector of 7 aero states for each element.
% Tar             : Transform from airfoil to rotorplane coordinates.
%                   3-by-3*Nel.
% Tas             : Transform from airfoil to section coordinates.
%                   3-by-3*Nel.
% TsB             : Transform from section to body coordinates.
%                   3-by-3*Nel.
% TBB0            : Transform from deformed to undeformed body coordinates.
%                   3-by-3*Nel, because we don't assume upfront the
%                   association between element and body.
% TB0g            : Transform from undeformed body to global coordinates.
%                   3-by-3*Nel.
% dTar            : Derivatives of transform from airfoil to rotorplane.
%                   3-by-3*24*Nel.  qy,qB,qn1,qn2.
% dTsB            : Derivatives of setion-to-body transform.
%                   3-by-3*12*Nel.  qn1,qn2.
% Vg              : 3*Nel, incoming wind, global coordinates.
% wg              : 3*Nel, element structural velocity, global coordinates.
% zr              : 2*Nel vector, x,y projections of element positions onto
%                   the rotorplane.
% ch,Lel          : Chord and lengths of the elements.  Length is measured 
%                   along the blade, that is, the distance between the 
%                   nodes.
% aoas,kfoils     : Splined airfoil coefficient tables.
% foilwt          : Nfoil-by-Nel table.  Weights to use when computing
%                   airfoil coefficients from splined tables.
% aoaz            : Zero-lift angles-of-attack for each element. Should
%                   be computed precisely from the airfoil tables.
% aoast           : 2*Nel vector, containing deep-stall angles-of-attack
%                   for each element.  Alternating positive, negative.
% xas,yas         : X^a and Y^a coordinates of the reference section
%                   coordinate system.
% rho             : Air density.
% Dp              : Diameter of the projected outer node on the rotorplane.
% rp              : Projected radius on the rotorplane.
% Lp              : Projected element length, not to be confused with Lel.
% Omega           : Rotor speed.
% 
%
% Outputs:
% --------
% A,B,C,D         : State matrices.

Nx = size(x,1);
Nel = Nx/7;
Ny = 137*Nel + 5;

Nb = 3;

if (psiFlag == 1)
   Neb = Nel/Nb;
end

% Allocate in full here, for speed, then sparsify at the end.
dxdt = zeros (Nx,1);
A    = zeros (Nx,Nx);
By   = zeros (Nx,Ny);
C    = zeros (Ny,Nx);
Dy   = zeros (Ny,Ny);

% Store the quantities that will be needed in the dynamic induced velocity
% calculation.
iViha  = zeros (3,1);
iVztsa = zeros (3,1);
iViya  = zeros (3,1);
ira    = zeros (3,1);
Vztsa  = zeros (9,1);
Viya   = zeros (6,1);
xera   = zeros (9,1);
xa     = zeros (12,1);
ra     = zeros (3,1);

% Reference indices (one before) for the final y variables.
iazi = Ny - 5;
iW   = Ny - 4;
iDp0 = Ny - 3;

elc = 0;
ann = 1;
for kel = 1:Nel

   elc = elc + 1;

   if (psiFlag == 0)
      iel = kel;
      iDp = iDp0;
   elseif (psiFlag == 1)
      iel = Neb*(elc-1) + ann;  % Index the Nb elements about an annulus.
      iDp = iDp0 + elc - 1;
   end

   ic72  = 72*(iel-1);
   ic36  = 36*(iel-1);
   ic9   =  9*(iel-1);
   ic3   =  3*(iel-1);
   ic2   =  2*(iel-1);

   [iad,ia1,ia2,iVih,iVi,                     ...
    iq,idq,ixng1,ixng2,ivng1,ivng2,           ...
    iwg,iwa,iaq,iC,iFldm,iFa,iFp,iFr,iFzts,   ...
    iVg,iUa,iUmag,iUr,iUzts,iVzts,iViq,       ...
    iViy,iVixyz,iWmag,ixeg,ixhg,              ...
    ixnr1,ixnr2,ixer,ir,iL,iz,iPr] = getis (iel);

   xer = [zr(ic2+[1:2]);0];  % Some functions below were written with xer as a 3-element
                             % vector.  The third entry isn't used.
   psie = atan2c (zr(ic2+2),zr(ic2+1));
   cpe = cos(psie);
   spe = sin(psie);

   % Transform the global structural motion to airfoil coordinates.
   [wa,ddy,Tag,dTag] = globalToAirfoil (wg(ic3+[1:3]),                        ...
                                        Tas(:,ic3+[1:3]),TsB(:,ic3+[1:3]),    ...
                                        TBB0(:,ic3+[1:3]),TB0g(:,ic3+[1:3]),  ...
                                        dTsB(:,ic36+[1:36]),dTBB0(:,ic9+[1:9]));
   Dy(iwa+[1:3],iq+[7:24]) = ddy(:,1:18);
   Dy(iwa+[1:3],iwg+[1:3]) = ddy(:,19:21);

   % Relate the xyz and zts components of induced velocity.
   Viz = x(iVi+1);
   Vit = x(iVi+2);
   Vis = 0;          % Let Vir_s = 0.  This should match BEMNL.
   Vizts = [Viz;Vit;Vis];
   [psie,Vixyz,ddy] = ztsToxyz (xer,Vizts);
   Dy(iVixyz+[1:3],ixer+[1:3]) = ddy(:,1:3); 
    C(iVixyz+[1:3],iVi+[1:2])  = ddy(:,4:5);

   % Compute the local flow velocity in airfoil coordinates.
   [Ua,ddu,ddy] = localFlow (Tag,dTag,Tar(:,ic3+[1:3]),dTar(:,ic72+[1:72]), ...
                             Vg(ic3+[1:3]),Vixyz,wa);
   Dy(iUa+[1:3],iq+[1:24])    = ddy(:,1:24);
   Dy(iUa+[1:3],iwa+[1:3])    = ddy(:,25:27);
   Dy(iUa+[1:3],iVixyz+[1:3]) = ddy(:,28:30);
   Dy(iUa+[1:3],iVg+[1:3])    = ddu;

   % Magnitude of local flow velocity.
   [Umag,ddy] = Umagnitude (Ua);
   Dy(iUmag+1,iUa+[1:2]) = ddy;

   % Local flow in the rotorplane CS.
   [Ur,ddy] = localFlowRotorCS (Tar(:,ic3+[1:3]),dTar(:,ic72+[1:72]),Ua);
   Dy(iUr+[1:3],iq+[1:24]) = ddy(:,1:24);
   Dy(iUr+[1:3],iUa+[1:3]) = ddy(:,25:27);

   % zts components of the local flow at the rotorplane.
   [psie,Uzts,ddy] = rotorUZTS (xer,Ur);
   Dy(iUzts+[1:3],ixer+[1:3]) = ddy(:,1:3);
   Dy(iUzts+[1:3],iUr+[1:3])  = ddy(:,4:6);

   % QS angle-of-attack.
   [aq,ddy] = alphaq (Ua);
   Dy(iaq+1,iUa+[1:3]) = ddy;

   % Dynamic angle-of-attack.
   [dxdt(iad+[1:3]),aa,bby] = aoadyn (ch(iel),x(iad+[1:3]),aq,aoast(ic2+[1:2]),Umag);
   A(iad+[1:3],iad+[1:3]) = aa;
   By(iad+[1:3],iaq+1)    = bby(:,1);
   By(iad+[1:3],iUmag+1)  = bby(:,2);

   % Airfoil coefficients.
   Caq = airfoilCoefficientsSpline (aoas,kfoils,foilwt(:,iel),aq);
   Cad = airfoilCoefficientsSpline (aoas,kfoils,foilwt(:,iel),x(iad+1));
   Co = Caq;
   Co(1) = Cad(1);
   Co(4) = Cad(4);
   [Cld,cc,ddy] = ClCdCm (Co,x(iad+1),aq,aoaz(iel));
   C(iC+[1:3],iad+1)  = cc;
   Dy(iC+[1:3],iaq+1) = ddy;

   % Airfoil forces, in airfoil coordinates.
   [Fld,Fa,ddy] = airfoilForces (rho,ch(iel),Lel(iel),xas(iel),yas(iel), ...
                                 aq,Umag,Cld,Co(2),Co(3));
   Dy(iFldm+[1:3],iaq+1)    = ddy(1:3,1);
   Dy(iFldm+[1:3],iUmag+1)  = ddy(1:3,2);
   Dy(iFldm+[1:3],iC+[1:3]) = ddy(1:3,3:5);
   Dy(iFa+[1:6],iaq+1)      = ddy(4:9,1);
   Dy(iFa+[1:6],iUmag+1)    = ddy(4:9,2);
   Dy(iFa+[1:6],iC+[1:3])   = ddy(4:9,3:5);

   % Airfoil forces in pitch coordinates.
   [Fp,ddy] = airfoilToPitch (Fa,TsB(:,ic3+[1:3]), ...
                              dTsB(:,ic36+[1:36]),Tas(:,ic3+[1:3]));
   Dy(iFp+[1:6],iq+[13:24]) = ddy(:,1:12);
   Dy(iFp+[1:6],iFa+[1:6])  = ddy(:,13:18);

   % zts components of the local flow at the rotorplane, not including
   % blade motion.
   [Vzts,ddy] = rotorVZTS (Uzts,x(iVi+[1:2]),rp(iel),Omega);
   Dy(iVzts+[1:3],iUzts+[1:3]) = ddy(:,1:3);
    C(iVzts+[1:3],iVi+[1:2])   = ddy(:,4:5);
   Dy(iVzts+[1:3],ir+1)        = ddy(:,6);
   Dy(iVzts+[1:3],iW+1)        = ddy(:,7);

   % Prandtl factor.
   [Pr,z,ddy] = prandtlLin (Nb,rp(iel),Dp(iel),Uzts);
   Dy(iz+1,iUzts+[1:3]) = ddy(2,[1:3]);
   Dy(iz+1,ir+1)        = ddy(2,4);
   Dy(iz+1,iDp+1)       = ddy(2,5);
   Dy(iPr+1,iz+1)       = ddy(1,6);

   % Transforming airfoil forces to the rotor CS.
   [Fr,ddy] = forcesRotorCS (Fa,Tar(:,ic3+[1:3]),dTar(:,ic72+[1:72]));
   Dy(iFr+[1:6],iq+[1:24]) = ddy(:,1:24);
   Dy(iFr+[1:6],iFa+[1:6]) = ddy(:,25:30);

   % Forces in zts coordinates, rotorUZTS works for this too.
   [psie,Fzts,ddy] = rotorUZTS (xer,Fr);
   Dy(iFzts+[1:3],ixer+[1:3]) = ddy(:,1:3);
   Dy(iFzts+[1:3],iFr+[1:3])  = ddy(:,4:6);

   % Quasi-steady induced velocity.
   [Viq,cc,ddy] = Vi_qs (3,x(iVi+1),Vzts,rho,rp(iel),Lp(iel),Pr,Fzts);
    C(iViq+[1:2],iVi+1)       = cc;
   Dy(iViq+[1:2],iVzts+[1:3]) = ddy(:,1:3);
   Dy(iViq+[1:2],ir+1)        = ddy(:,4);
   Dy(iViq+[1:2],iL+1)        = ddy(:,5);
   Dy(iViq+[1:2],iPr+1)       = ddy(:,6);
   Dy(iViq+[1:2],iFzts+[1:3]) = ddy(:,7:9);

   % Rotorplane velocity magnitude.
   [Wmag,cc,ddy] = Wvel (x(iVi+1),Vzts,Pr);
    C(iWmag+1,iVi+1)       = cc;
   Dy(iWmag+1,iVzts+[1:3]) = ddy(:,1:3);
   Dy(iWmag+1,iPr+1)       = ddy(:,4);

   % Vi corrected for yaw angle.
   [Viy,cc,ddy] = Viyaw (x(iVi+1),Viq,Vzts,rp(iel),Dp(iel),Pr,Wmag);
    C(iViy+[1:2],iVi+1)       = cc;
   Dy(iViy+[1:2],iVzts+[1:3]) = ddy(:,1:3);
   Dy(iViy+[1:2],iViq+[1:2])  = ddy(:,4:5);
   Dy(iViy+[1:2],ir+1)        = ddy(:,6);
   Dy(iViy+[1:2],iDp+1)       = ddy(:,7);
   Dy(iViy+[1:2],iPr+1)       = ddy(:,8);
   Dy(iViy+[1:2],iWmag+1)     = ddy(:,9);

   if (psiFlag == 1)
      % Store annulus variables needed to call Vidyn.
      e4 = 4*(elc-1);
      e3 = 3*(elc-1);
      e2 = 2*(elc-1);
      iViha(elc)      = iVih;
      iVztsa(elc)     = iVzts;
      iViya(elc)      = iViy;
      ira(elc)        = ir;
      Vztsa(e3+[1:3]) = Vzts;
      Viya(e2+[1:2])  = Viy;
      xera(e3+[1:3])  = xer;
      xa(e4+[1:4])    = x(iVih+[1:4]);
      ra(elc)         = rp(iel);
   end

   if (psiFlag == 0)
      [dxdt(iVih+[1:4]),aa,bby] = Vidyn (psiFlag,x(iVih+[1:4]),Vzts,Viy, ...
                                         rp(iel),Dp(iel),azi,Omega);
       A(iVih+[1:4],iVih+[1:4])  = aa;
      By(iVih+[1:4],iVzts+[1:3]) = bby(:,1:3);
      By(iVih+[1:4],iViy+[1:2])  = bby(:,4:5);
      By(iVih+[1:4],ir+1)        = bby(:,6);
      By(iVih+[1:4],iazi+1)      = bby(:,7);
      By(iVih+[1:4],iW+1)        = bby(:,8);
      By(iVih+[1:4],iDp+1)       = bby(:,9);
   elseif ((psiFlag == 1) && (elc == 3))
      % Finished with an annulus, now call Vidyn and fill matrices.
      % [Vidyn verified via complex step.]
      elc = 0;
      ann = ann + 1;
      [dx,aa,bby] = Vidyn (psiFlag,xa,Vztsa,Viya,ra,Dp,azi,Omega);
      for jj = 1:3
         jc4 = 4*(jj-1);
          A(iViha(jj)+[1:4],iViha(1)+[1:4])  =  aa(jc4+[1:4],1:4);
          A(iViha(jj)+[1:4],iViha(2)+[1:4])  =  aa(jc4+[1:4],5:8);
          A(iViha(jj)+[1:4],iViha(3)+[1:4])  =  aa(jc4+[1:4],9:12);
         By(iViha(jj)+[1:4],iVztsa(1)+[1:3]) = bby(jc4+[1:4],1:3);
         By(iViha(jj)+[1:4],iViya(1)+[1:2])  = bby(jc4+[1:4],4:5);
         By(iViha(jj)+[1:4],ira(1)+1)        = bby(jc4+[1:4],6);
         By(iViha(jj)+[1:4],iVztsa(2)+[1:3]) = bby(jc4+[1:4],7:9);
         By(iViha(jj)+[1:4],iViya(2)+[1:2])  = bby(jc4+[1:4],10:11);
         By(iViha(jj)+[1:4],ira(2)+1)        = bby(jc4+[1:4],12);
         By(iViha(jj)+[1:4],iVztsa(3)+[1:3]) = bby(jc4+[1:4],13:15);
         By(iViha(jj)+[1:4],iViya(3)+[1:2])  = bby(jc4+[1:4],16:17);
         By(iViha(jj)+[1:4],ira(3)+1)        = bby(jc4+[1:4],18);
         By(iViha(jj)+[1:4],iazi+1)          = bby(jc4+[1:4],19);
         By(iViha(jj)+[1:4],iW+1)            = bby(jc4+[1:4],20);
         By(iViha(jj)+[1:4],iDp0+1)          = bby(jc4+[1:4],21);
         By(iViha(jj)+[1:4],iDp0+2)          = bby(jc4+[1:4],22);
         By(iViha(jj)+[1:4],iDp0+3)          = bby(jc4+[1:4],23);

         dxdt(iViha(jj)+[1:4]) = dx(jc4+[1:4]);
      end

   end

end

A  = sparse(A);
By = sparse(By);
C  = sparse(C);
Dy = sparse(Dy);

