function [LHS,A,Bu,By,C,Du,Dy] =                                    ...
                  linkAeroelastic (s,a,                             ...
                                   lls,aas,bbus,bbys,ccs,ddus,ddys, ...
                                   lla,aaa,bbya,cca,ddya,           ...
                                   ddypv,ddyRSA)
%
% In many studies, the systems and controls may vary in architecture,
% requiring non-standard linking procedures; while the turbine
% aeroelastic architecture remains consistent.  It is thus convenient
% to have a standard function that links the aerodynamic and 
% structural parts.
%
%   States:              y vector:             u vector:
% ----------------------- Structure --------------------------
%   eta      Neta        q         Ndj         F          Ndj
%   deta/dt  Neta        dq/dt     Ndj         d2eta/dt2  Neta
%                        d2q/dt2   Ndj
%                        xng     3*Nnod
%                        vng     6*Nnod
%                        F         Ndj
% ---------------------- Aerodynamic -------------------------
%   ad         1         (q)       1:24        Vg          3
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
%                        r         124
%                        Lp        135
%                        z         136
%                        f         137  
% (Repeat the above x,y,u consecutively for each blade element.)
%                        Azi       (1)
%                        Waero     (1)  (aero rotor speed)
%                        Dp        (1)
%
% Version:        Changes:
% --------        -------------
% 12.04.2018      Original code.
% 22.01.2020      Minor modifications to account for d2eta/dt2 as
%                 part of the u vector.
%
% Version:        Verification:
% --------        -------------
% 12.04.2018      
% 22.01.2020      
%
% Inputs:
% -------
% 
%
% Outputs:
% --------
% A ... Dy        : State matrices.

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
Ndj = Ndof + 6;
Nnod = Ndof/6;
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);

Nxs  = size(aas,1);
Nus  = size(bbus,2);
Nys1 = size(ccs,1);
Nys  = Nys1 + 9*Nnod + Ndj;
N    = Nxs/2;

Nxa  = size(aaa,1);
Nya  = size(cca,1);
Nua  = 3*a.Nb*a.Neb;  % We will define Vg as a global input.
Nyae = 137;           % Number of aero y's associated with each element.
Nea  = a.Nb*a.Neb;
Neb  = a.Neb;

% ====================================================
% Fill the state space with each module
Nx    = Nxs + Nxa;
Ny    = Nys + Nya;
Nu    = Nus + Nua;
nnzL  = nnz(lls) + nnz(lla);
nnzA  = round(1.2*(nnz(aas)  + nnz(aaa)));
nnzBu = round(2.0*(nnz(bbus)));
nnzBy = round(1.2*(nnz(bbys) + nnz(bbya) + nnz(bbus)));
nnzC  = round(1.2*(nnz(ccs)  + nnz(cca)));
nnzDu = round(2.0*(nnz(ddus)));
nnzDy = round(1.2*(nnz(ddys) + nnz(ddya) + nnz(ddus)));
LHS   = spalloc (Nx,Nx,nnzL);
A     = spalloc (Nx,Nx,nnzA);
Bu    = spalloc (Nx,Nu,nnzBu);
By    = spalloc (Nx,Ny,nnzBy);
C     = spalloc (Ny,Nx,nnzC);
Du    = spalloc (Ny,Nu,nnzDu);
Dy    = spalloc (Ny,Ny,nnzDy);

LHS(1:Nxs,1:Nxs)   = lls;
  A(1:Nxs,1:Nxs)   = aas;
 Bu(1:Nxs,1:Nus)   = bbus;
 By(1:Nxs,1:Nys1)  = bbys;
 By(1:Nxs,Nys1+9*Nnod+[1:Ndj]) = bbus(:,1:Ndj);  % Aero forces in y vector.
  C(1:Nys1,1:Nxs)  = ccs;
 Du(1:Nys1,1:Nus)  = ddus;
 Dy(1:Nys1,1:Nys1) = ddys;
 Dy(1:Nys1,Nys1+9*Nnod+[1:Ndj]) = ddus(:,1:Ndj); % Aero forces in y vector.

LHS(Nxs+[1:Nxa],Nxs+[1:Nxa]) = lla;
  A(Nxs+[1:Nxa],Nxs+[1:Nxa]) = aaa;
 By(Nxs+[1:Nxa],Nys+[1:Nya]) = bbya;
  C(Nys+[1:Nya],Nxs+[1:Nxa]) = cca;
 Dy(Nys+[1:Nya],Nys+[1:Nya]) = ddya;

% ====================================================
% Links between modules

% Get aero element DOF references.
[iad,ia1,ia2,iVih,iVi,                   ...
 iq,idq,ixng1,ixng2,ivng1,ivng2,         ...
 iwg,iwa,iaq,iC,iFldm,iFa,iFp,iFr,iFzts, ...
 iVg,iUa,iUmag,iUr,iUzts,iVzts,iViq,     ...
 iViy,iVixyz,iWmag,ixeg,ixhg,            ...
 ixnr1,ixnr2,ixer,irp,iLp,izPr,iPr] = getis (1);

% Link the aero element qy, qB, qn1, and qn2, as well as the 
% corresponding velocities, to the structural values.
for icomp = 1:6

   indr = Nys + iq + icomp + [0:Nyae:Nyae*(Nea-1)].';
   indc = idofs(3) + icomp;
   Dy(indr,indc) = 1;

   indr = Nys + iq + 6 + icomp + [0:Nyae:Nyae*(Nea-1)].';
   indc = icomp + [idofs(6) idofs(7) idofs(8)].';
   mat = sparse(Nea,3);
   mat(1:Neb,1) = 1;
   mat(Neb+[1:Neb],2) = 1;
   mat(2*Neb+[1:Neb],3) = 1;
   Dy(indr,indc) = mat;

   indr = Nys + iq + 12 + icomp + [0:Nyae:Nyae*(Nea-1)].';
   indc = icomp + [idofs(6)+[0:6:6*(Neb-1)] ...
                   idofs(7)+[0:6:6*(Neb-1)] ...
                   idofs(8)+[0:6:6*(Neb-1)]].';
%   Dy(indr,indc) = speye(Nea);  % Should not link qn1 to qB at first node.
   Dy(indr,indc) = [sparse(1,Nea);                                      ...
                    sparse(Neb-1,1) speye(Neb-1) sparse(Neb-1,2*Neb);   ...
                    sparse(1,Nea);                                      ...
                    sparse(Neb-1,Neb+1) speye(Neb-1) sparse(Neb-1,Neb); ...
                    sparse(1,Nea);                                      ...
                    sparse(Neb-1,2*Neb+1) speye(Neb-1)];

   indr = Nys + iq + 18 + icomp + [0:Nyae:Nyae*(Nea-1)].';
   indc = icomp + 6 + [idofs(6)+[0:6:6*(Neb-1)] ...
                       idofs(7)+[0:6:6*(Neb-1)] ...
                       idofs(8)+[0:6:6*(Neb-1)]].';
   Dy(indr,indc) = speye(Nea);

   indr = Nys + idq + icomp + [0:Nyae:Nyae*(Nea-1)].';
   indc = Ndj + idofs(3) + icomp;
   Dy(indr,indc) = 1;

   indr = Nys + idq + 6 + icomp + [0:Nyae:Nyae*(Nea-1)].';
   indc = Ndj + icomp + [idofs(6) idofs(7) idofs(8)].';
   mat = sparse(Nea,3);
   mat(1:Neb,1) = 1;
   mat(Neb+[1:Neb],2) = 1;
   mat(2*Neb+[1:Neb],3) = 1;
   Dy(indr,indc) = mat;

   indr = Nys + idq + 12 + icomp + [0:Nyae:Nyae*(Nea-1)].';
   indc = Ndj + icomp + [idofs(6)+[0:6:6*(Neb-1)] ...
                         idofs(7)+[0:6:6*(Neb-1)] ...
                         idofs(8)+[0:6:6*(Neb-1)]].';
%   Dy(indr,indc) = speye(Nea);
   Dy(indr,indc) = [sparse(1,Nea);                                      ...
                    sparse(Neb-1,1) speye(Neb-1) sparse(Neb-1,2*Neb);   ...
                    sparse(1,Nea);                                      ...
                    sparse(Neb-1,Neb+1) speye(Neb-1) sparse(Neb-1,Neb); ...
                    sparse(1,Nea);                                      ...
                    sparse(Neb-1,2*Neb+1) speye(Neb-1)];

   indr = Nys + idq + 18 + icomp + [0:Nyae:Nyae*(Nea-1)].';
   indc = Ndj + icomp + 6 + [idofs(6)+[0:6:6*(Neb-1)] ...
                             idofs(7)+[0:6:6*(Neb-1)] ...
                             idofs(8)+[0:6:6*(Neb-1)]].';
   Dy(indr,indc) = speye(Nea);

end

% Link xng and vng, which are needed in the aero calculations.
indr = Nys1 + [1:9*Nnod].';
indc = [[1:Ndof] Ndj+[1:Ndof]].';
Dy(indr,indc) = Dy(indr,indc) + ddypv;

% Link the aero element xhg, xng, vng, and wg to the structural values.
for icomp = 1:3

   indr = Nys + ixhg  + icomp + [0:Nyae:Nyae*(Nea-1)].';
   indc = 3*Ndj + 3*(inodm(6) - 2) + icomp;
   Dy(indr,indc) = 1;

   indr = Nys + ixng1 + icomp + [0:Nyae:Nyae*(Nea-1)].';
   indc = 3*Ndj + icomp                   ...
        + [3*(inods(6)-1)+[0:3:3*(Neb-1)] ...
           3*(inods(7)-1)+[0:3:3*(Neb-1)] ...
           3*(inods(8)-1)+[0:3:3*(Neb-1)]].';
   Dy(indr,indc) = speye(Nea);

   indr = Nys + ixng2 + icomp + [0:Nyae:Nyae*(Nea-1)].';
   indc = 3*Ndj + 3 + icomp               ...
        + [3*(inods(6)-1)+[0:3:3*(Neb-1)] ...
           3*(inods(7)-1)+[0:3:3*(Neb-1)] ...
           3*(inods(8)-1)+[0:3:3*(Neb-1)]].';
   Dy(indr,indc) = speye(Nea);

   indr = Nys + ivng1 + icomp + [0:Nyae:Nyae*(Nea-1)].';
   indc = 3*Ndj + 3*Nnod + icomp          ...
        + [6*(inods(6)-1)+[0:6:6*(Neb-1)] ...
           6*(inods(7)-1)+[0:6:6*(Neb-1)] ...
           6*(inods(8)-1)+[0:6:6*(Neb-1)]].';
   Dy(indr,indc) = speye(Nea);

   indr = Nys + ivng2 + icomp + [0:Nyae:Nyae*(Nea-1)].';
   indc = 3*Ndj + 3*Nnod + 6 + icomp      ...
        + [6*(inods(6)-1)+[0:6:6*(Neb-1)] ...
           6*(inods(7)-1)+[0:6:6*(Neb-1)] ...
           6*(inods(8)-1)+[0:6:6*(Neb-1)]].';
   Dy(indr,indc) = speye(Nea);

end

% Link the aerodynamic rotor speed.  Let the aerodynamic rotor speed
% be defined as the z^d rotational speed of the shaft at the hub
% center.  Note that the control measurement, in the control section
% below, is the difference between the shaft speed and nacelle motion.
% The two are therefore not exactly equal.
indr  = Nys + Nyae*Nea + [1:2].';
indc  = [1:2*Ndj].';
Dy(indr,indc) = Dy(indr,indc) + ddyRSA;

% Make Vg a global input.
for icomp = 1:3

   ind3 = Nus       + icomp + [0:3:3*(Nea-1)].';
   indc = Nys + iVg + icomp + [0:Nyae:Nyae*(Nea-1)].';
   Bu(:,ind3) = By(:,indc);
   Du(:,ind3) = Dy(:,indc);

end

% Link aero forces to structural forces.  Note that the airfoil forces
% are given for the blade elements, while the structural forces are
% applied to the nodes.
for icomp = 1:6

   for ib = 1:3

      indr = 3*Ndj + 9*Nnod + idofs(5+ib) + icomp + [0:6:6*(Neb-1)].';
      indc = Nys + Nyae*Neb*(ib-1) + iFp + icomp + [0:Nyae:Nyae*(Neb-1)].';
      Dy(indr,indc) = 0.5*speye(Neb);

      indr = 3*Ndj + 9*Nnod + idofs(5+ib) + icomp + [6:6:6*Neb].';
      Dy(indr,indc) = Dy(indr,indc) + 0.5*speye(Neb);

   end

end

% Link input forces to output forces, so that these represent the
% full force vector.  Zero the direct link between input forces
% and the equations of motion.  When the matrices are linked
% later, this will nonetheless produce the correct B matrix.
indr = 3*Ndj + 9*Nnod+ [1:Ndj].';
indc = [1:Ndj].';
Du(indr,indc) = speye(Ndj);
Bu(1:Nxs,indc) = sparse(Nxs,Ndj);
